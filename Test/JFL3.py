import random
import numpy as np
import multiprocessing
import time

import Mapper

#Глобальная переменная 
MAP_PARAM = None


def process_bed_files(file1, file2):
    def process_lines(lines):
        return [(int(parts[1]), int(parts[2])) for parts in 
                (line.split() for line in lines if line.startswith("chr")) 
                if len(parts) >= 3 and parts[1].isdigit() and parts[2].isdigit()]
    
    def process_lines2(lines):
        return [float(parts[3]) for parts in 
                (line.split() for line in lines if line.startswith("chr")) 
                if len(parts) > 3 and parts[3].replace('.', '', 1).isdigit()]
    
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1, lines2 = f1.readlines(), f2.readlines()

    return process_lines(lines1), process_lines(lines2), process_lines2(lines2)


def map_genes_to_domains(Domens, Genes, expressions):
    if MAP_PARAM == 1:
        # 1. Отнести ген целиком к тому домену, в котором находится большая часть гена
        arr_of_domens = Mapper.map_genes_to_domains_1(Domens, Genes, expressions)
        
    elif MAP_PARAM == 2:
        # 2. Отнести экспрессию гена пропорционально длине, попадающей в домен
        arr_of_domens = Mapper.map_genes_to_domains_2(Domens, Genes, expressions)
    
    elif MAP_PARAM == 3:
        # 3. Отнести экспрессию гена целиком к точке старта
        arr_of_domens = Mapper.map_genes_to_domains_3(Domens, Genes, expressions)
    
    elif MAP_PARAM == 4:
        # 4. Отнести экспрессию гена целиком к точке конца транскрипции
        arr_of_domens = Mapper.map_genes_to_domains_4(Domens, Genes, expressions)
    else: raise ValueError("Некорректный режим в map_param.txt. Ожидалось 1, 2, 3 или 4.")

    return arr_of_domens

    
def rotate_genes(Domens, Genes, expressions, mean_gene_len):
    shift_value = np.random.randint(0, mean_gene_len) * np.random.choice([-1, 1])
    chrom_length = max(gene[1] for gene in Genes)

    new_expressions = [None] * len(Genes)
    new_genes = []
    
    for i, (start, end) in enumerate(Genes):
        new_start = (start + shift_value) % chrom_length
        new_end = (end + shift_value) % chrom_length
        
        if new_end < new_start:
            mode = MAP_PARAM
            if mode == 1 or mode == 2:
                if new_end > chrom_length - new_start:
                    new_start = 0
                else:
                    new_end = chrom_length
            elif mode == 3:
                new_end = chrom_length
            elif mode == 4:
                new_start = 0
            else:
                raise ValueError(f'Ошибка поворота - мод мапера указан неверно: {mode}')
        
        new_genes.append((new_start, new_end))
        new_expressions[i] = expressions[i]
    arr_of_domens = map_genes_to_domains(Domens, new_genes, new_expressions)
    return arr_of_domens

def z_criterium(mu1: float, mu2: float, sig1: float, sig2: float, n: int) -> float:
    return (mu1 - mu2) / np.sqrt(sig1/n + sig2/n)

def pad_and_stack(arr_of_domens):
    max_length = max(len(domen) if domen is not None else 0 for domen in arr_of_domens)
    padded_array = np.full((len(arr_of_domens), max_length), np.nan)
    for i, domen in enumerate(arr_of_domens):
        if domen is not None:
            padded_array[i, :len(domen)] = domen
    return padded_array

def count_statvalues_horis(arr_of_domens: np.ndarray) -> tuple:
    arr_of_domens = arr_of_domens[~np.isnan(arr_of_domens).all(axis=1)]
    overall_mean = np.nanmean(arr_of_domens)
    overall_dispersion = np.nanvar(arr_of_domens, ddof=1) 
    # Чтобы дисперсия рассчитывалась по формуле с поправкой на выборку
    return overall_mean, overall_dispersion

def one_chromosome_rework_new(Domens, Genes, expressions, iterations) -> list:
    global MAP_PARAM
    if MAP_PARAM is None:
        with open("map_param.txt", "r") as f:
            MAP_PARAM = int(f.read().strip())
    
    mean_gene_len = sum(end - start for start, end in Genes) / len(Genes) if Genes else 0
    mu1, d1 = count_statvalues_horis(pad_and_stack(map_genes_to_domains(Domens, Genes, expressions)))
    
    mu_list, d_list, matrix, mu_v, d_vertical = [], [], [], [], []
    
    for _ in range(iterations):
        new_arr_of_domens = rotate_genes(Domens, Genes, expressions, mean_gene_len)
        matrix.append(new_arr_of_domens)
        
    for new_arr_of_domens in matrix:
        mu2, d2 = count_statvalues_horis(pad_and_stack(new_arr_of_domens))
        if mu2 is not None:
            mu_list.append(mu2)
        if d2 is not None:
            d_list.append(d2)
            
    mu_new = np.mean(mu_list) if mu_list else 0
    d_new = np.mean(d_list) if d_list else 0

    transposed_matrix = np.transpose(matrix)    
    for new_arr_of_domens in transposed_matrix:
        filtered_arr = [x for x in new_arr_of_domens if x is not None]
        if len(filtered_arr) > 0:
            mu2, d = count_statvalues_horis(pad_and_stack(filtered_arr))
            #print(mu2, d2)
            if mu2 is not None:
                mu_v.append(mu2)
            if d is not None:
                d_vertical.append(d)
    #print(d_vertical)        
    mu_vert = np.mean(mu_v) if mu_v else 0
    d_vert = np.mean(d_vertical)
    z = z_criterium(mu1, mu_new, d1, d_new, len(Genes)) - 1
    result = [z, mu1, d1, mu_new, d_new, mu_vert, d_vert]
    return result
