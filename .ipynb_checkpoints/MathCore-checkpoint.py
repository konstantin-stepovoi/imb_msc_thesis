import numpy as np
import DomenSorter
from parameters import Config

'''
этот модуль в программе основной. он абстрагирован от того, какие данные мы ему скармливаем
он принимает уже готовый arr_of_domens и отдаёт z_criterium
arr_of_domens может поступать как от parser, так и от modeler+domensorter
'''
config = Config()

def count_statvalues(arr_of_domens: np.ndarray) -> tuple:
    means = []
    variances = []

    for domen in arr_of_domens:
        if domen is not None and len(domen) > 0:
            means.append(np.mean(domen))
            variances.append(np.var(domen))
    
    overall_mean = np.mean(means) if means else None
    overall_dispersion = np.mean(variances) if variances else None

    return overall_mean, overall_dispersion

def z_criterium(mu1: float, mu2: float, sig1: float, sig2: float, n: int) -> float:
    return (mu1 - mu2) / (sig1 + sig2) * np.sqrt(n)


def one_chromosome_rework(chromosome) -> float:
    '''
    вот тут надо добавить флажок с какими данными мы работаем - принятыми из Modeler
    или из Parser.
    но это понадобится делать только тогда, когда parser будет готов
    '''

    arr_of_domens = DomenSorter.Sort_to_domens(chromosome, config.domen_length)
    mu1, d1 = count_statvalues(arr_of_domens)
    mu_list, d_list = [], []
    for _ in range(config.rotations_per_sample):
        shift_value = int(np.random.rand() * config.mean_gene_length)
        chromosome.rotate(shift_value)
        arr_of_domens = DomenSorter.Sort_to_domens(chromosome, config.domen_length)
        mu2, d2 = count_statvalues(arr_of_domens)
        if mu2 is not None:
            mu_list.append(mu2)
        if d2 is not None:
            d_list.append(d2) 
    mu2 = np.mean(mu_list)
    d2 = np.mean(d_list)
    return z_criterium(mu1, mu2, d1, d2, config.number_of_genes)