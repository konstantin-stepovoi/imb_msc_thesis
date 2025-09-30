#!/usr/bin/python3
import numpy as np
import multiprocessing
import time
import logging
from numba import jit
import Domain_mapper
from Domain_mapper import *

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def get_mapper_method():
    global _mapper_method_cache
    try:
        return _mapper_method_cache
    except NameError:
        pass

    filename = 'cur_metadata.txt'
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.strip().startswith('mapper_method'):
                    parts = line.strip().split(' ', maxsplit=1)
                    if len(parts) == 2:
                        value = parts[1].strip().lower()
                        if value not in {'proportional', 'start', 'stop', 'middle'}:
                            raise ValueError(f"Unsupported mapper method: {value}")
                        _mapper_method_cache = value
                        return _mapper_method_cache
                    else:
                        raise ValueError("Incorrect syntax in metadata: expected 'mapper_method -- value'")
    except FileNotFoundError:
        raise RuntimeError(f"Файл {filename} не найден.")
    raise RuntimeError("Parameter 'mapper_method' not found in cur_metadata.txt.")



def process_bed_files(file1, file2, metadata_file='cur_metadata.txt'):
    def read_analysis_type(path):
        with open(path) as f:
            for line in f:
                parts = line.strip().split(":", 1)
                if len(parts) == 2 and parts[0].strip().lower() == "analysis_type":
                    return parts[1].strip().lower()
        return 'chromosome-wide'  # fallback default

    def parse_regions(lines, shift_needed):
        regions = []
        shift = 0
        prev_start = -1
        for line in lines:
            if not line.startswith("chr"):
                continue
            parts = line.split()
            if len(parts) < 3 or not parts[1].isdigit() or not parts[2].isdigit():
                continue
            start, end = int(parts[1]), int(parts[2])
            if shift_needed and start < prev_start:
                shift += prev_start  # новая хромосома, сдвигаем
            regions.append((start + shift, end + shift))
            prev_start = start
        return regions

    def parse_scores(lines, shift_needed):
        scores = []
        shift = 0
        prev_start = -1
        for line in lines:
            if not line.startswith("chr"):
                continue
            parts = line.split()
            if len(parts) >= 4 and parts[1].isdigit():
                start = int(parts[1])
                if shift_needed and start < prev_start:
                    shift += prev_start
                try:
                    scores.append(float(parts[3]))
                except ValueError:
                    continue
                prev_start = start
        return scores

    analysis_type = read_analysis_type(metadata_file)
    shift_needed = (analysis_type == 'genome-wide')

    with open(file1) as f1, open(file2) as f2:
        lines1, lines2 = f1.readlines(), f2.readlines()

    return (
        parse_regions(lines1, shift_needed),
        parse_regions(lines2, shift_needed),
        parse_scores(lines2, shift_needed)
    )


##############################################################################################################
# РАЗНЫЕ СПОСОБЫ МАПИРОВАНИЯ. МОЖЕТ ПОТОМ ВЫБРОСИТЬ В ОТДЕЛЬНЫЙ МОДУЛЬ?
# ------------------------------
# Global, lightweight state for rotations
# ------------------------------

class _State:
    def __init__(self):
        self.mapper_method: str = ""
        self.domain_lengths: List[int] = []
        self.flat_expr: Optional[np.ndarray] = None  # concatenated once
        self.total_count: int = 0
        # Proportional-specific
        self.dom_weights: Optional[List[np.ndarray]] = None  # per-domain weights
        self.flat_expr_overlaps: Optional[np.ndarray] = None  # concatenated raw expr per-overlap
        self.starts: Optional[np.ndarray] = None

STATE = _State()

# ------------------------------
# Unified mapper with state init
# ------------------------------

def map_genes_to_domains(domains, genes, expressions, mapper_method: str):
    """
    Returns list[np.ndarray] per domain and initializes STATE for rotations.
    mapper_method in {"start","stop","middle","proportional"}
    """
    exprs = Domain_mapper._ensure_exprs(genes, expressions)
    STATE.mapper_method = mapper_method

    if mapper_method == "start":
        dom_arrays = map_domains_start(domains, genes, exprs)
    elif mapper_method == "stop":
        dom_arrays = map_domains_stop(domains, genes, exprs)
    elif mapper_method == "middle":
        dom_arrays = map_domains_center(domains, genes, exprs)
    elif mapper_method == "proportional":
        idxs, ws = build_overlap_weight_lists(domains, genes)
        # Convert weights to ndarrays once (saves overhead later)
        ws_nd = [np.asarray(w, dtype=float) for w in ws]
        dom_arrays = map_domains_weighted_sparse(exprs, idxs, ws)

        # Build flat per-overlap expression vector aligned to weights order
        flat_idx = np.fromiter((gi for sub in idxs for gi in sub), dtype=int)
        flat_expr_overlaps = exprs[flat_idx] if flat_idx.size else np.empty(0, dtype=float)

        # Persist proportional-state
        STATE.dom_weights = [w.copy() for w in ws_nd]
        STATE.flat_expr_overlaps = flat_expr_overlaps
    else:
        raise ValueError(f"Unknown mapper method: {mapper_method}")

    # Persist non-proportional flat view (concatenate ONCE)
    STATE.domain_lengths = [len(a) for a in dom_arrays]
    STATE.starts = np.cumsum([0] + STATE.domain_lengths[:-1]).astype(int)
    if mapper_method != "proportional":
        STATE.flat_expr = np.concatenate(dom_arrays) if STATE.starts.size and sum(STATE.domain_lengths) else np.empty(0, dtype=float)
        STATE.total_count = STATE.flat_expr.size
    else:
        STATE.flat_expr = None
        STATE.total_count = STATE.flat_expr_overlaps.size

    return dom_arrays


##############################################################################################################
#Основная математика

def _fast_mean_var(arr: np.ndarray):
    if arr.size == 0:
        return 0.0, 0.0
    m = float(arr.sum() / arr.size)
    v = float(((arr - m) ** 2).sum() / arr.size)
    return m, v

def count_statvalues_horis_no_pad(domains_arrays: list[np.ndarray]):
    """
    Рассчитывает среднее и дисперсию по доменам, **фильтруя пустые и одиночные домены**.

    Аргументы:
    ----------
    domains_arrays : list[np.ndarray]
        Список массивов экспрессий по доменам. Каждый элемент соответствует одному домену.

    Возвращает:
    ---------
    overall_mean : float
        Среднее значение по полноценным доменам (с >=2 генами).
    overall_dispersion : float
        Средняя дисперсия по полноценным доменам.
    """
    # Фильтруем пустые и домены с 1 геном
    filtered_domains = [d for d in domains_arrays if len(d) > 1]

    if not filtered_domains:
        return 0.0, 0.0

    # Средние и дисперсии по каждому домену
    means = np.array([d.mean() for d in filtered_domains], dtype=float)
    vars_ = np.array([d.var() for d in filtered_domains], dtype=float)

    # Итоговые показатели
    overall_mean = float(means.mean())
    overall_dispersion = float(vars_.mean())

    return overall_mean, overall_dispersion

'''
СТАРАЯ ВЕРСИЯ
def count_statvalues_horis_no_pad(domains_arrays: List[np.ndarray]):
    means = []
    vars_ = []
    for d in domains_arrays:
        m, v = _fast_mean_var(d)
        means.append(m)
        vars_.append(v)
    means = np.asarray(means, dtype=float)
    vars_ = np.asarray(vars_, dtype=float)
    overall_mean = float(means.mean()) if means.size else 0.0
    overall_dispersion = float(vars_.mean()) if vars_.size else 0.0
    return overall_mean, overall_dispersion'''


def z_criterium(mu1: float, mu2: float, sig1: float, sig2: float, n: int) -> float:
    denom = sig1 / max(n,1) + sig2 / max(n,1)
    if denom <= 0:
        return 0.0
    return float((mu1 - mu2) / np.sqrt(denom))

def _lazy_slice_cyclic(base: np.ndarray, start: int, length: int) -> np.ndarray:
    """Return a contiguous view of length `length` starting at `start` with modulo wrap.
    This avoids np.roll on the whole array; it creates at most one concatenation per slice
    only when wrapping is needed.
    """
    n = base.size
    if n == 0 or length == 0:
        return np.empty(0, dtype=base.dtype)
    s = start % n
    e = s + length
    if e <= n:
        return base[s:e]
    # wrap-around: stitch two views (small, local copy only for this domain)
    first = base[s:]
    second = base[:e - n]
    return np.concatenate((first, second))


def rotate_genes(current_shift: int) -> List[np.ndarray]:
    """
    Produce per-domain arrays for the given cyclic shift WITHOUT rebuilding global arrays.
    - For non-proportional: returns slices of the once-concatenated flat_expr.
    - For proportional: fetch raw expr per-overlap then multiply by precomputed weights per domain.
    Memory: O(sum(domain_lengths)) temporaries per call, no quadratic blow-up, no global growth.
    """
    if STATE.total_count == 0:
        return [np.empty(0, dtype=float) for _ in STATE.domain_lengths]

    out: List[np.ndarray] = []
    starts = STATE.starts
    lens = STATE.domain_lengths

    if STATE.mapper_method == "proportional":
        base = STATE.flat_expr_overlaps
        assert base is not None and STATE.dom_weights is not None
        pos = current_shift % base.size
        for L, w in zip(lens, STATE.dom_weights):
            raw = _lazy_slice_cyclic(base, pos, L)
            out.append(raw * w)
            pos += L
        return out

    # Non-proportional branches
    base = STATE.flat_expr
    pos = current_shift % base.size
    for L in lens:
        out.append(_lazy_slice_cyclic(base, pos, L))
        pos += L
    return out


# ------------------------------
# Main routine (drop-in replacement)
# ------------------------------

def Weilford_function(Domens, Genes, expressions, iterations, mapper_method: str) -> list:
    if mapper_method not in {"middle", "proportional", "start", "stop"}:
        raise ValueError(f"unknown mapping method: {mapper_method}")

    # Init mapping + state
    original_domains = map_genes_to_domains(Domens, Genes, expressions, mapper_method)

    # Baseline metrics
    mu1, d1 = count_statvalues_horis_no_pad(original_domains)

    # Iterations (no global growth)
    total_sum = 0.0
    total_sumsq = 0.0
    total_count = 0

    mu_list_sum = 0.0
    d_list_sum = 0.0
    valid_iter_count = 0

    rng = np.random.default_rng()

    for i in range(iterations):
        # random shift in [-N, N)
        shift = int(rng.integers(-STATE.total_count, STATE.total_count)) if STATE.total_count else 0
        rotated = rotate_genes(shift)

        if (i + 1) % 10 == 0:
            logging.info(f"Iteration {i+1}/{iterations} completed.")

        mu2, d2 = count_statvalues_horis_no_pad(rotated)
        mu_list_sum += mu2
        d_list_sum += d2
        valid_iter_count += 1

        # Accumulate population stats
        for dom in rotated:
            if dom.size == 0:
                continue
            total_sum += float(dom.sum())
            total_sumsq += float((dom * dom).sum())
            total_count += int(dom.size)

    mu_new = mu_list_sum / max(valid_iter_count, 1)
    d_new = d_list_sum / max(valid_iter_count, 1)
    mu_final = total_sum / max(total_count, 1)
    d_final = (total_sumsq / max(total_count, 1)) - mu_final ** 2 if total_count else 0.0

    z = z_criterium(mu1, mu_new, mu_final, d_final, max(iterations, 1))

    mu_vert = 'Ok'
    d_vert = 'Ok'

    return [z, mu1, d1, mu_new, d_new, mu_vert, d_vert, mu_final, d_final]



# ------------------------------
# Convenience wrapper
# ------------------------------

def run_with_config(domains, genes, expressions, iterations, mapper_method: str):
    return Weilford_function(domains, genes, expressions, iterations, mapper_method)