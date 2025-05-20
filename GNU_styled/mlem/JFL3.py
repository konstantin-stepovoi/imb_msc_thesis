#!/usr/bin/python3
import random
import numpy as np
import multiprocessing
import time

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



def process_bed_files(file1, file2):
    def parse_regions(lines):
        return [
            (int(parts[1]), int(parts[2]))
            for line in lines if line.startswith("chr")
            if (parts := line.split()) and len(parts) >= 3
            and parts[1].isdigit() and parts[2].isdigit()
        ]

    def parse_scores(lines):
        scores = []
        for line in lines:
            if line.startswith("chr"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        scores.append(float(parts[3]))
                    except ValueError:
                        continue
        return scores

    with open(file1) as f1, open(file2) as f2:
        lines1, lines2 = f1.readlines(), f2.readlines()

    return parse_regions(lines1), parse_regions(lines2), parse_scores(lines2)

##############################################################################################################
# РАЗНЫЕ СПОСОБЫ МАПИРОВАНИЯ. МОЖЕТ ПОТОМ ВЫБРОСИТЬ В ОТДЕЛЬНЫЙ МОДУЛЬ?
        
def map_genes_to_domains_start(Domens, Genes, expressions): 
    arr_of_domens = np.empty(len(Domens), dtype=object)
    for domain_index, (dom_start, dom_end) in enumerate(Domens):
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            if dom_start <= gene_start < dom_end:
                this_domen.append(expressions[i])
        arr_of_domens[domain_index] = np.array(this_domen)
    return arr_of_domens

def map_genes_to_domains_middle(Domens, Genes, expressions): 
    arr_of_domens = np.empty(len(Domens), dtype=object)
    for domain_index, (dom_start, dom_end) in enumerate(Domens):
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            gene_middle = (gene_start + gene_end) // 2
            if dom_start <= gene_middle < dom_end:
                this_domen.append(expressions[i])
        arr_of_domens[domain_index] = np.array(this_domen)
    return arr_of_domens

def map_genes_to_domains_stop(Domens, Genes, expressions): 
    arr_of_domens = np.empty(len(Domens), dtype=object)
    for domain_index, (dom_start, dom_end) in enumerate(Domens):
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            if dom_start < gene_end <= dom_end:
                this_domen.append(expressions[i])
        arr_of_domens[domain_index] = np.array(this_domen)
    return arr_of_domens

def map_genes_to_domains_propor(Domens, Genes, expressions): 
    arr_of_domens = np.empty(len(Domens), dtype=object)
    gene_index = 0
    for domain_index, (dom_start, dom_end) in enumerate(Domens):
        this_domen = []
        while gene_index < len(Genes) and Genes[gene_index][1] < dom_start:
            gene_index += 1
        temp_index = gene_index  
        while temp_index < len(Genes) and Genes[temp_index][0] < dom_end:
            gene_start, gene_end = Genes[temp_index]
            gene_length = gene_end - gene_start
            expression = expressions[temp_index]
            if dom_start <= gene_start and gene_end <= dom_end:
                this_domen.append(expression)
            elif gene_start <= dom_start and dom_end <= gene_end:
                overlap_length = dom_end - dom_start
                this_domen.append(expression * (overlap_length / gene_length))
            elif gene_start < dom_start < gene_end <= dom_end:
                overlap_length = gene_end - dom_start
                this_domen.append(expression * (overlap_length / gene_length))
            elif dom_start <= gene_start < dom_end < gene_end:
                overlap_length = dom_end - gene_start
                this_domen.append(expression * (overlap_length / gene_length))
            temp_index += 1
        arr_of_domens[domain_index] = np.array(this_domen)
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

def pad_and_stack(arr_of_domens):
    max_length = max(len(domen) for domen in arr_of_domens)
    padded_array = np.array([np.pad(domen, (0, max_length - len(domen)), constant_values=np.nan)
                             for domen in arr_of_domens])
    return padded_array
    
def count_statvalues_horis(arr_of_domens: np.ndarray) -> tuple:
    means = np.nanmean(arr_of_domens, axis=1)
    variances = np.nanvar(arr_of_domens, axis=1)

    overall_mean = np.nanmean(means)
    overall_dispersion = np.nanmean(variances)
    return overall_mean, overall_dispersion

def z_criterium(mu1: float, mu2: float, sig1: float, sig2: float, n: int) -> float:
    return (mu1 - mu2) / np.sqrt(sig1/n + sig2/n)

def rotate_genes(Domens, Genes, expressions, shift_value):
    # Находим общую длину хромосомы (конец последнего гена)
    chrom_length = max(gene[1] for gene in Genes)

    # Сдвигаем координаты всех генов
    new_genes = []
    for start, end in Genes:
        new_start = (start + shift_value) % chrom_length
        new_end = (end + shift_value) % chrom_length
        new_genes.append((new_start, new_end))

    # Массив expressions должен быть также сдвинут
    new_expressions = expressions[-shift_value % len(expressions):] + expressions[:-shift_value % len(expressions)]

    # Вызываем функцию маппинга с новыми сдвинутыми данными
    arr_of_domens = map_genes_to_domains(Domens, new_genes, new_expressions)
    filtered_domens = pad_and_stack(arr_of_domens)
    return filtered_domens


def one_chromosome_rework_new(Domens, Genes, expressions, iterations) -> list:
    mapper_method = get_mapper_method()
    if mapper_method not in {'middle', 'proportional', 'start', 'stop'}:
        raise ValueError(f'unknown type of mapping method in md file: {mapper_method}')
    
    gene_lengths = [end - start for start, end in Genes]
    mean_gene_len = sum(gene_lengths) / len(gene_lengths) if gene_lengths else 0

    # Исходные горизонтальные метрики
    mu1, d1 = count_statvalues_horis(pad_and_stack(map_genes_to_domains(Domens, Genes, expressions)))

    mu_list_sum = 0
    d_list_sum = 0
    valid_iter_count = 0

    # Вертикальное накопление
    num_domens = len(Domens)
    vertical_sums = [0.0 for _ in range(num_domens)]
    vertical_sumsq = [0.0 for _ in range(num_domens)]
    vertical_counts = [0 for _ in range(num_domens)]

    for _ in range(iterations):
        shift_value = int(np.random.rand() * mean_gene_len)
        shift_value *= np.random.choice([-1, 1])

        new_arr_of_domens = rotate_genes(Domens, Genes, expressions, shift_value)

        # Горизонтальные метрики
        mu2, d2 = count_statvalues_horis(new_arr_of_domens)
        if mu2 is not None:
            mu_list_sum += mu2
            d_list_sum += d2
            valid_iter_count += 1

        # Вертикальное накопление
        for j, dom_expr in enumerate(new_arr_of_domens):
            dom_expr = np.array(dom_expr)
            dom_expr = dom_expr[~np.isnan(dom_expr)]
            if len(dom_expr) == 0:
                continue
            vertical_sums[j] += np.sum(dom_expr)
            vertical_sumsq[j] += np.sum(dom_expr ** 2)
            vertical_counts[j] += len(dom_expr)

    # Финальное горизонтальное
    mu_new = mu_list_sum / valid_iter_count if valid_iter_count else 0
    d_new = d_list_sum / valid_iter_count if valid_iter_count else 0

    # Z критерий
    z = z_criterium(mu1, mu_new, d1, d_new, len(Genes)) - 1

    # Финальные вертикальные метрики
    vertical_means = []
    vertical_vars = []

    for j in range(num_domens):
        n = vertical_counts[j]
        if n == 0:
            continue
        mu = vertical_sums[j] / n
        var = (vertical_sumsq[j] / n) - mu**2
        vertical_means.append(mu)
        vertical_vars.append(var)

    mu_vert = np.mean(vertical_means) if vertical_means else 0
    d_vert = np.mean(vertical_vars) if vertical_vars else 0

    return [z, mu1, d1, mu_new, d_new, mu_vert, d_vert]
