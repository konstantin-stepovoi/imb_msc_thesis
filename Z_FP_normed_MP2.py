#!/usr/bin/python3
#Импорты
import random
from dataclasses import dataclass, field
import numpy as np
from typing import Optional, List
from scipy.integrate import quad
from scipy.interpolate import interp1d
import concurrent.futures
import multiprocessing
import time

# Определяем параметры модели
y0 = 2.18266
x0 = -356.69431
A1 = 15019.02824
t1 = 93.14008
A2 = 16085.18644
t2 = 93.13557
A3 = 67.24575
t3 = 2165.97117

# Определяем диапазон генерации случайных чисел
xmin = 50
xmax = 66000

# Функция тройного экспоненциального затухания
def exp_decay3(x, y0, x0, A1, t1, A2, t2, A3, t3):
    return (
        y0
        + A1 * np.exp(-(x - x0) / t1)
        + A2 * np.exp(-(x - x0) / t2)
        + A3 * np.exp(-(x - x0) / t3)
    )

# Нормализация функции (интегрируемая плотность вероятности)
integral, _ = quad(exp_decay3, xmin, xmax, args=(y0, x0, A1, t1, A2, t2, A3, t3))

def normalized_exp_decay3(x, y0, x0, A1, t1, A2, t2, A3, t3):
    return exp_decay3(x, y0, x0, A1, t1, A2, t2, A3, t3) / integral

# Вычисляем CDF путем интегрирования нормализованной PDF
x_values = np.linspace(xmin, xmax, 1000)
cdf_values = np.cumsum(normalized_exp_decay3(x_values, y0, x0, A1, t1, A2, t2, A3, t3))
cdf_values /= cdf_values[-1]  # нормализуем CDF к 1

# Исправляем CDF для предотвращения ошибки при интерполяции
inverse_cdf = interp1d(cdf_values, x_values, bounds_error=False, fill_value=(xmin, xmax))

def generate_random_exp_decay3(n=1):
    random_uniform_values = np.random.rand(n)  # генерация n случайных значений от 0 до 1
    random_uniform_values = np.clip(random_uniform_values, cdf_values.min(), cdf_values.max())
    return inverse_cdf(random_uniform_values)

#основные датаклассы для последующего хранения данных:
@dataclass
class Gene:
    length: np.int32  #randint(860, 190001)
    expression: np.uint8 #np.random.randint(0, 50)
    start_point: np.int32
    end_point: np.int32
    name: Optional[str] = field(default = None)


    def __post_init__(self):
        # Дополнительная валидация или вычисления могут быть добавлены здесь.
        if self.length != self.end_point - self.start_point:
            raise ValueError("Length does not match the difference between start and end points.")

class Chromosome(np.ndarray):
    def __new__(cls, genes: List[Gene] = None):
        if genes is None:
            genes = []
        obj = np.asarray(genes, dtype=object).view(cls)
        return obj
    
    def add_gene(self, gene: Gene):
        new_size = len(self) + 1
        new_array = np.empty(new_size, dtype=object)  # Создаем новый массив нужного размера
        new_array[:-1] = self  # Копируем данные из старого массива
        new_array[-1] = gene  # Добавляем новый ген
        return new_array.view(Chromosome)
    
    def index(self, search_value: int, by: str = 'start_point'):
        for i, gene in enumerate(self):
            if hasattr(gene, by) and getattr(gene, by) == search_value:
                return i
        raise ValueError(f"{search_value} not found in Chromosome by {by}")

    def rotate(self, shift_value: int):
        chromosome_length = self[-1].end_point
        for gene in self:
            gene.start_point += shift_value
            if gene.start_point > chromosome_length:
                gene.start_point = gene.start_point - chromosome_length
            gene.end_point += shift_value
            if gene.end_point > chromosome_length:
                gene.end_point = gene.end_point - chromosome_length

        indices = np.argsort([gene.start_point for gene in self])
        self[:] = self[indices]
        return self
def Chromosome_generator(num: int):
    chromosome = Chromosome()
    coordinate = 0
    for j in range(num):
        len = int(generate_random_exp_decay3(1)[0])
        gene = Gene(length = len, expression=np.random.randint(0, 50), 
                    start_point = coordinate, end_point = coordinate + len)
        coordinate = gene.end_point
        chromosome = chromosome.add_gene(gene)
    #print(coordinate)
        
    return chromosome

def create_empty_object_array(size):
    return np.empty(size, dtype=object)

def add_array(array_of_arrays, index, new_array):
    array_of_arrays[index] = new_array

def Sort_to_domens(chromosome, limit):
    # Определите количество доменов
    number_of_domens = (chromosome[-1].end_point // limit) + 1  # Плюс один для учёта последнего домена
    arr_of_domens = create_empty_object_array(number_of_domens)
    
    for n in range(number_of_domens):
        this_domen = []
        this_d_l = limit * n
        this_d_r = limit * (n + 1)
    
        for gene in chromosome:
            gene_start_in_domen = gene.start_point % limit
            gene_end_in_domen = gene.end_point % limit
            
            # Если ген полностью в текущем домене
            if this_d_l <= gene.start_point < this_d_r and this_d_l < gene.end_point <= this_d_r:
                this_domen.append(gene.expression)
            # Если ген пересекает начало домена
            elif gene.start_point < this_d_l < gene.end_point:
                overlap_length = gene.end_point - this_d_l
                this_domen.append(gene.expression * overlap_length / gene.length)
            # Если ген пересекает конец домена
            elif this_d_l <= gene.start_point < this_d_r and gene.end_point > this_d_r:
                overlap_length = this_d_r - gene.start_point
                this_domen.append(gene.expression * overlap_length / gene.length)
            # Если ген полностью вне домена (как слева, так и справа)
            elif gene.end_point < this_d_l or gene.start_point >= this_d_r:
                continue  # Пропускаем гены, которые не попадают в текущий домен
    
        if len(this_domen) > 0:
            add_array(arr_of_domens, n, np.array(this_domen))
    
    return arr_of_domens


def count_statvalues(arr_of_domens):
    means = []
    variances = []

# Шаг 1: Для каждого подмассива считаем среднее и дисперсию
    for domen in arr_of_domens:
        if domen is not None and len(domen) > 0:  # Проверяем, что подмассив не пустой
            mean = np.mean(domen)
            variance = np.var(domen)
            means.append(mean)
            variances.append(variance)
            
    if means and variances: 
        overall_mean = np.mean(means)
        overall_dispersion = np.mean(variances)
    else:
        overall_mean = None
        overall_dispersion = None

    return(overall_mean, overall_dispersion)


def z_criterium(mu1, mu2, sig1, sig2, n):
    return (mu1 - mu2)/(sig1 + sig2) * n**0.5

#от тут надо будет вынести потом  limit  в глобал, это мед длин дом

def one_chromosome_rework():
    limit = 110000
    mu = []
    d = []
    chromosome = Chromosome_generator(5000)
    arr_of_domens = Sort_to_domens(chromosome, limit)
    mu1, d1 = count_statvalues(arr_of_domens)
    for i in range(10000):
        shift_value = np.random.rand()*34000
        chromosome.rotate(int(shift_value))
        arr_of_domens = Sort_to_domens(chromosome, limit)
        m2, d2 = count_statvalues(arr_of_domens)
        if m2 is not None: mu.append(m2)
        if d2 is not None: d.append(d2)
    mu2 = np.mean(mu)
    #print(len(mu))
    d2 = np.mean(d)
    #print(len(d))
    return z_criterium(mu1, mu2, d1, d2, 5000)

def append_number_to_file(number, file_path):
    try:
        with open(file_path, 'a') as file:
            file.write(str(number) + '\n')  # Добавляем число и перевод строки
        #print(f"Число {number} успешно добавлено в файл {file_path}.")
    except Exception as e:
        print(f"Произошла ошибка при записи в файл: {e}")

file_path = str(input('filepath:'))

def worker_func(j):
    # Используем тайм стэмп для генерации уникального сида
    seed = int((time.time() * 1e6) % 1e6) + j
    np.random.seed(seed)
    random.seed(seed)
    
    z = one_chromosome_rework()
    print(f'try {j}, z = {z}')
    return z  # Отправляем результат в очередь

def writer_func(results, file_path):
	with open(file_path, 'a') as file:
		for z in results:
			file.write(str(z) + '\n')  # Пишем результат в файл

def run_in_processes_with_writer(num_tasks):
    # Используем Pool для запуска рабочих процессов
    with multiprocessing.Pool() as pool:
        # Передаем задачи в пул процессов
        results = pool.map(worker_func, range(num_tasks))

    # Когда все процессы завершены, отправляем сигнал завершения записи
    writer_process = multiprocessing.Process(target = writer_func, args=(results, file_path))
    writer_process.start()
    writer_process.join()  # Ждём завершения процесса записи

# Запуск задач в многопроцессорном режиме
if __name__ == '__main__':
    run_in_processes_with_writer(int(input('Start N:')))
