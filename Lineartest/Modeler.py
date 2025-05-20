from Gene_containers import Gene, Chromosome
from scipy.integrate import quad
from scipy.interpolate import interp1d
from parameters import Config
from random import randint
import numpy as np
from datetime import datetime

'''
этот модуль только создаёт новые хромосомы, ориентируясь на
Gene_containers и на parameters. 
потом от отдаёт хромосому в DomenSorter на сортировку генов по доменам
'''
config = Config()

# Функция тройного экспоненциального затухания
def exp_decay3(x, y0, x0, A1, t1, A2, t2, A3, t3) -> float:
    return (
        y0
        + A1 * np.exp(-(x - x0) / t1)
        + A2 * np.exp(-(x - x0) / t2)
        + A3 * np.exp(-(x - x0) / t3)
    )

# Нормализация функции
integral, _ = quad(
    exp_decay3,
    config.Min_gene_len,
    config.Max_gene_len,
    args=(config.y0, config.x0, config.A1, config.t1, config.A2, config.t2, config.A3, config.t3)
)

# Нормализованная функция затухания
def normalized_exp_decay3(x: np.ndarray) -> np.ndarray:
    return exp_decay3(
        x, config.y0, config.x0, config.A1, config.t1, config.A2, config.t2, config.A3, config.t3
    ) / integral

# Вычисление CDF
x_values = np.linspace(config.Min_gene_len, config.Max_gene_len, 1000)
cdf_values = np.cumsum(normalized_exp_decay3(x_values))
cdf_values /= cdf_values[-1]

# Создание обратной CDF
inverse_cdf = interp1d(cdf_values, x_values, 
                       bounds_error=False, fill_value=(config.Min_gene_len, config.Max_gene_len))

# Генерация случайных значений на основе тройного экспоненциального затухания
def generate_random_exp_decay3(n: int = 1) -> np.ndarray:
    random_uniform_values = np.random.rand(n)
    random_uniform_values = np.clip(random_uniform_values, cdf_values.min(), cdf_values.max())
    return inverse_cdf(random_uniform_values)

def uniform_distribution() -> int:
    random_l = randint(config.Min_gene_len, config. Max_gene_len)
    return random_l

def Chromosome_generator(num: int) -> Chromosome:
    chromosome = Chromosome()
    coordinate = 0
    for _ in range(num):
        gene_length = int(generate_random_exp_decay3(1)[0])
        gene = Gene(
            length=gene_length,
            expression=np.random.randint(config.minimal_expression, config.maximal_expression),
            start_point=coordinate,
            end_point=coordinate + gene_length
        )
        coordinate = gene.end_point
        chromosome = chromosome.add_gene(gene)
    return chromosome

'''
форматирование должно работать вот так:
chr1    1000    2000    Gene1   .   +
chr1    3000    4500    Gene2   .   -
chr1    5000    6000    Gene3   .   +
chr1    6500    7000    Gene4   .   -
chr1    8000    9000    Gene5   .   +
'''

def Generate_Chromosome_f(number_chroms: int, chrom_length: int):
    filename = 'generated_chromosomes_' + datetime.now().strftime('%y_%m_%d_%H_%M_%S') + '.txt'
    
    with open(filename, 'a') as file:
        for chroms in range(number_chroms):
            array_lengths = []
            current_chrom_length = chrom_length  # Сохраняем длину для каждой хромосомы
            while current_chrom_length > 0:
                if config.genetype == 'uni':
                    gene_l = uniform_distribution()
                elif config.genetype == 'exp':
                    gene_l = int(generate_random_exp_decay3(1)[0])
                else:
                    raise ValueError("Invalid type: choose 'uni' or 'exp'")
                    
                array_lengths.append(gene_l)
                current_chrom_length -= gene_l
            
            current_position = 0
            for current_gene in range(1, len(array_lengths) + 1):
                expression = np.random.randint(config.minimal_expression, config.maximal_expression)
                file.write(f'chr{chroms} {current_position} {current_position + array_lengths[current_gene - 1]} Gene{current_gene} {expression} +\n')
                current_position += array_lengths[current_gene - 1]
    return filename
        