#!/usr/bin/python3
import numpy as np
import multiprocessing
import time
import logging
from numba import jit
import Mapping_types as Mp

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
                if line.strip().startswith('Seed'):
                    seed = line.strip().split()[1]
                    np.random.seed(seed)
            
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


CALIBRATION_WEIGHTS = [] 
GENE_COUNT = 0
        
def map_genes_to_domains_start(Domens, Genes, expressions): 
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    for dom_start, dom_end in Domens:
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            if dom_start <= gene_start < dom_end:
                this_domen.append(expressions[i])

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.ones(len(this_domen), dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens


def map_genes_to_domains_middle(Domens, Genes, expressions): 
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    for dom_start, dom_end in Domens:
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            gene_middle = (gene_start + gene_end) // 2
            if dom_start <= gene_middle < dom_end:
                this_domen.append(expressions[i])

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.ones(len(this_domen), dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens

def map_genes_to_domains_stop(Domens, Genes, expressions): 
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    for dom_start, dom_end in Domens:
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            if dom_start < gene_end <= dom_end:
                this_domen.append(expressions[i])

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.ones(len(this_domen), dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens

    
def map_genes_to_domains_propor(Domens, Genes, expressions):
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    gene_index = 0
    for dom_start, dom_end in Domens:
        this_domen = []
        this_weights = []

        while gene_index < len(Genes) and Genes[gene_index][1] < dom_start:
            gene_index += 1
        temp_index = gene_index

        while temp_index < len(Genes) and Genes[temp_index][0] < dom_end:
            gene_start, gene_end = Genes[temp_index]
            gene_length = gene_end - gene_start
            expression = expressions[temp_index]

            weight = 0.0
            if dom_start <= gene_start and gene_end <= dom_end:
                weight = 1.0
            elif gene_start <= dom_start and dom_end <= gene_end:
                weight = (dom_end - dom_start) / gene_length
            elif gene_start < dom_start < gene_end <= dom_end:
                weight = (gene_end - dom_start) / gene_length
            elif dom_start <= gene_start < dom_end < gene_end:
                weight = (dom_end - gene_start) / gene_length

            if weight > 0:
                this_domen.append(expression * weight)
                this_weights.append(weight)

            temp_index += 1

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.array(this_weights, dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens


def map_genes_to_domains(Domens, Genes, expressions):
    mapper_method = get_mapper_method()
    if mapper_method == "proportional":
        return map_genes_to_domains_propor(Domens, Genes, expressions)
    elif mapper_method == "start":
        return map_genes_to_domains_start(Domens, Genes, expressions)
    elif mapper_method == "stop":
        return map_genes_to_domains_stop(Domens, Genes, expressions)
    elif mapper_method == "middle":
        return map_genes_to_domains_middle(Domens, Genes, expressions)
    else:
        raise ValueError(f"Unknown mapper method: {mapper_method}")


def count_statvalues_horis_no_pad(arr_of_domens: np.ndarray) -> tuple:
    def fast_mean(arr):
        return arr.sum() / arr.size

    def fast_var(arr):
        mu = fast_mean(arr)
        return ((arr - mu) ** 2).sum() / arr.size

    means = np.array([fast_mean(domen) for domen in arr_of_domens])
    variances = np.array([fast_var(domen) for domen in arr_of_domens])

    overall_mean = means.mean()
    overall_dispersion = variances.mean()

    return overall_mean, overall_dispersion

def z_criterium(mu1: float, mu2: float, sig1: float, sig2: float, n: int) -> float:
    if sig2 == 0:
        return 0.0 if sig1 == 0 else 10.0

    # Эмпирическая компенсация: при больших n немного корректируем отношение дисперсий
    correction_factor = 1 + 0.2 * (1 - np.exp(-n / 1000))  # При n=0 → 1, при n→∞ → 1.2
    adj_sig1 = sig1
    adj_sig2 = sig2 * correction_factor
    Q = Mp.affine_mod(sig1, sig2)
    Z = max(adj_sig1, adj_sig2) / min(adj_sig1, adj_sig2)
    direction = 1 if adj_sig1 > adj_sig2 else -1
    score = Q #direction * np.log(Z)
    return score


def rotate_genes(arr_of_domens):
    global CALIBRATION_WEIGHTS, GENE_COUNT

    flat_expr = np.concatenate(arr_of_domens)
    shift = int(np.random.uniform(-1, 1) * GENE_COUNT)
    rotated_expr = np.roll(flat_expr, shift)

    new_arr = np.empty(len(CALIBRATION_WEIGHTS), dtype=object)
    index = 0
    for i, weights in enumerate(CALIBRATION_WEIGHTS):
        n = len(weights)
        new_arr[i] = rotated_expr[index: index + n] * weights
        index += n
    return new_arr


def Weilford_old(Domens, Genes, expressions, iterations) -> list:
    mapper_method = get_mapper_method()
    if mapper_method not in {'middle', 'proportional', 'start', 'stop'}:
        raise ValueError(f'unknown type of mapping method in md file: {mapper_method}')


    # Исходные домены
    original_domains = map_genes_to_domains_propor(Domens, Genes, expressions)
    lengths = np.array([len(sublist) for sublist in original_domains])

    mu1, d1 = count_statvalues_horis_no_pad(original_domains)

    # Инициализация накоплений
    mu_list_sum = 0
    d_list_sum = 0
    valid_iter_count = 0

    new_arr_of_domens = original_domains
    for i in range(iterations):
        new_arr_of_domens = rotate_genes(new_arr_of_domens)
        mu2, d2 = count_statvalues_horis_no_pad(new_arr_of_domens)
        if mu2 is not None:
            mu_list_sum += mu2
            d_list_sum += d2
            valid_iter_count += 1

    # Финальные метрики по итерациям
    mu_new = mu_list_sum / valid_iter_count if valid_iter_count else 0
    d_new = d_list_sum / valid_iter_count if valid_iter_count else 0

    # Откалиброванный z - test
    z = z_criterium(mu1, mu_new, d1, d_new, iterations)

    return [z, mu1, d1, mu_new, d_new, 'Ok', 'Ok']


def Weilford_function(Domens, Genes, expressions, iterations):
    mapper_method = get_mapper_method()
    if mapper_method not in {'middle', 'proportional', 'start', 'stop'}:
        raise ValueError(f'unknown type of mapping method in md file: {mapper_method}')

        
    orig = map_genes_to_domains(Domens, Genes, expressions)
    orig_domains = orig
    
    mean_before, var_before = count_statvalues_horis_no_pad(orig_domains)
    global CALIBRATION_WEIGHTS
    length = len(CALIBRATION_WEIGHTS)
    
    means = np.empty(iterations)
    vars = np.empty(iterations)
    
    rotated = rotate_genes(orig_domains)
    
    for num in range(iterations):
        rotated = rotate_genes(orig_domains)
        mean_new, var_new = count_statvalues_horis_no_pad(rotated)
        
        means[num] = mean_new
        vars[num] = var_new
    fin_mean = np.sum(means)/iterations
    fin_var = np.sum(vars)/iterations
    z = z_criterium(mean_before, fin_mean, var_before, fin_var, iterations)
    return [z, mean_before, fin_mean, var_before, fin_var, 'Ok', 'Ok']