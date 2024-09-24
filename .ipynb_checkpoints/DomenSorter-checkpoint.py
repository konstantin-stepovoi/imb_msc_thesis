import numpy as np
from Gene_containers import Chromosome

'''
этот файл принимает хромосому и отдаёт массив доменов,
каждый из элементов которого заполнен согласно определяемыми 
нами правилами и файлу параметров
'''
def create_empty_object_array(size: int) -> np.ndarray:
    return np.empty(size, dtype=object)

def add_array(array_of_arrays: np.ndarray, index: int, new_array: np.ndarray):
    array_of_arrays[index] = new_array

def Sort_to_domens(chromosome: Chromosome, limit: int) -> np.ndarray:
    number_of_domens = (chromosome[-1].end_point // limit) + 1
    arr_of_domens = create_empty_object_array(number_of_domens)
    
    for n in range(number_of_domens):
        this_domen = []
        this_d_l = limit * n
        this_d_r = limit * (n + 1)
    
        for gene in chromosome:
            gene_start_in_domen = gene.start_point % limit
            gene_end_in_domen = gene.end_point % limit
            
            if this_d_l <= gene.start_point < this_d_r and this_d_l < gene.end_point <= this_d_r:
                this_domen.append(gene.expression)
            elif gene.start_point < this_d_l < gene.end_point:
                overlap_length = gene.end_point - this_d_l
                this_domen.append(gene.expression * overlap_length / gene.length)
            elif this_d_l <= gene.start_point < this_d_r and gene.end_point > this_d_r:
                overlap_length = this_d_r - gene.start_point
                this_domen.append(gene.expression * overlap_length / gene.length)
    
        if this_domen:
            add_array(arr_of_domens, n, np.array(this_domen))
    
    return arr_of_domens
