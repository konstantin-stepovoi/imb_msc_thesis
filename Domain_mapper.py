from __future__ import annotations
from typing import List, Tuple, Optional
import numpy as np
import logging


"""
Модуль domain_mapper
--------------------

Набор функций для маппинга экспрессий генов на домены по разным критериям.

Единый интерфейс всех функций:
--------------------------------
map_domains_<method>(domains, genes, expressions=None, weights=None) -> List[List[float]]

Аргументы:
----------
domains : List[Tuple[int, int]]
    Список доменов [(start, stop), ...].

genes : List[Tuple[int, int, float]] или List[Tuple[int, int]]
    Список генов:
        - Если expressions=None: каждая запись — (start, stop, expression).
        - Если expressions передан: каждая запись — (start, stop).

expressions : List[float], optional
    Экспрессии генов в том же порядке, что и genes.
    Если None — берётся третье значение из кортежа genes.

weights : List[List[float]], optional
    Матрица весов [gene][domain] для метода weighted.
    Если None, будет рассчитана автоматически (см. build_overlap_weight_matrix).

Возвращает:
-----------
List[List[float]]
    Список списков экспрессий по доменам. Каждый подсписок соответствует одному домену.

Методы:
-------
map_domains_start(domains, genes, expressions=None)
    Маппинг по старту гена.

map_domains_stop(domains, genes, expressions=None)
    Маппинг по стопу гена.

map_domains_center(domains, genes, expressions=None)
    Маппинг по центру гена.

build_overlap_weight_matrix(domains, genes)
    Построение матрицы весов перекрытия генов и доменов.

map_domains_weighted(domains, genes, expressions=None, weights=None)
    Пропорциональный маппинг по перекрытию с весами.
"""


def _ensure_exprs(genes: List[Tuple[int, int, float]] | List[Tuple[int, int]], expressions: Optional[List[float]]):
    if expressions is not None:
        return np.asarray(expressions, dtype=float)
    # expect triplets (start, stop, expr)
    return np.asarray([g[2] for g in genes], dtype=float)


def map_domains_start(domains, genes, expressions=None) -> List[np.ndarray]:
    exprs = _ensure_exprs(genes, expressions)
    positions = np.asarray([(g[0], g[1]) for g in genes], dtype=int)
    out: List[np.ndarray] = []
    gene_idx = 0
    n = len(positions)
    for dom_start, dom_stop in domains:
        buf = []
        while gene_idx < n:
            g_start = positions[gene_idx, 0]
            if g_start > dom_stop:
                break
            if dom_start <= g_start <= dom_stop:
                buf.append(exprs[gene_idx])
            gene_idx += 1
        out.append(np.asarray(buf, dtype=float))
    return out


def map_domains_stop(domains, genes, expressions=None) -> List[np.ndarray]:
    exprs = _ensure_exprs(genes, expressions)
    positions = np.asarray([(g[0], g[1]) for g in genes], dtype=int)
    out: List[np.ndarray] = []
    gene_idx = 0
    n = len(positions)
    for dom_start, dom_stop in domains:
        buf = []
        while gene_idx < n:
            g_stop = positions[gene_idx, 1]
            if g_stop > dom_stop:
                break
            if dom_start <= g_stop <= dom_stop:
                buf.append(exprs[gene_idx])
            gene_idx += 1
        out.append(np.asarray(buf, dtype=float))
    return out


def map_domains_center(domains, genes, expressions=None) -> List[np.ndarray]:
    exprs = _ensure_exprs(genes, expressions)
    positions = np.asarray([(g[0], g[1]) for g in genes], dtype=int)
    out: List[np.ndarray] = []
    gene_idx = 0
    n = len(positions)
    for dom_start, dom_stop in domains:
        buf = []
        while gene_idx < n:
            g_start, g_stop = positions[gene_idx]
            g_center = (g_start + g_stop) // 2
            if g_center > dom_stop:
                break
            if dom_start <= g_center <= dom_stop:
                buf.append(exprs[gene_idx])
            gene_idx += 1
        out.append(np.asarray(buf, dtype=float))
    return out



def build_overlap_weight_matrix(domains, genes):
    weights = []
    for g_start, g_stop, *_ in genes:
        gene_len = g_stop - g_start
        row = []
        for dom_start, dom_stop in domains:
            overlap_start = max(g_start, dom_start)
            overlap_stop = min(g_stop, dom_stop)
            overlap_len = max(0, overlap_stop - overlap_start)
            weight = overlap_len / gene_len if gene_len > 0 else 0.0
            row.append(weight)
        weights.append(row)
    return weights


def map_domains_weighted(domains, genes, expressions=None, weights=None):
    exprs = expressions if expressions is not None else [g[2] for g in genes]
    if weights is None:
        weights = build_overlap_weight_matrix(domains, genes)

    num_domains = len(weights[0])
    domain_exprs = [[] for _ in range(num_domains)]

    for gene_idx, expr in enumerate(exprs):
        for dom_idx, w in enumerate(weights[gene_idx]):
            if w > 0:
                domain_exprs[dom_idx].append(w * expr)

    return domain_exprs


# ------------------------------
# Proportional mapping (sparse)
# ------------------------------

def build_overlap_weight_lists(domains: List[Tuple[int,int]],
                               genes: List[Tuple[int,int]] | List[Tuple[int,int,float]]):
    """
    Return per-domain sparse overlap structure:
      - dom_gene_idx: List[List[int]]  (indices of genes contributing to domain)
      - dom_weights:  List[List[float]] (corresponding fractional overlaps)
    Complexity: O(N_genes + N_domains) for mostly sorted inputs using two-pointer sweep.
    """
    dom_gene_idx: List[List[int]] = [[] for _ in range(len(domains))]
    dom_weights:  List[List[float]] = [[] for _ in range(len(domains))]

    # Ensure plain (start, stop) tuples for genes
    gene_coords = [(g[0], g[1]) for g in genes]

    # Two-pointer sweep over sorted intervals
    d = 0
    for gi, (gs, ge) in enumerate(gene_coords):
        if ge <= gs:
            continue
        # advance domain pointer until potential overlap
        while d < len(domains) and domains[d][1] <= gs:
            d += 1
        # walk overlapping domains
        j = d
        while j < len(domains):
            ds, de = domains[j]
            if ds >= ge:
                break
            # overlap
            overlap = max(0, min(ge, de) - max(gs, ds))
            if overlap > 0:
                dom_gene_idx[j].append(gi)
                dom_weights[j].append(overlap / (ge - gs))
            j += 1
    return dom_gene_idx, dom_weights


def map_domains_weighted_sparse(exprs: np.ndarray,
                                dom_gene_idx: List[List[int]],
                                dom_weights: List[List[float]]) -> List[np.ndarray]:
    out: List[np.ndarray] = []
    for idxs, ws in zip(dom_gene_idx, dom_weights):
        if idxs:
            e = exprs[idxs]
            w = np.asarray(ws, dtype=float)
            out.append(e * w)
        else:
            out.append(np.empty(0, dtype=float))
    return out
