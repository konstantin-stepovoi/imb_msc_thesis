#!/usr/bin/python3
# Импорты
# Interface.py
import numpy as np
import random, time, multiprocessing
from typing import List
from datetime import datetime
import Modeler
from parameters import Config
from MathCore import one_chromosome_rework
from Parser import load_real_data  # Предположим, что этот модуль будет позже написан
'''
это просто надстройка, которая спрашивает у нас, что надо сделать и 
пишет результаты в файл
'''
config = Config()
generated_data = ''

def append_number_to_file(number: float, file_path: str):
    try:
        with open(file_path, 'a') as file:
            file.write(f"{number}\n")
    except Exception as e:
        print(f"Error writing to file: {e}")

def write_parameters_to_file(file_path: str, config: Config):
    try:
        with open(file_path, 'w') as file:
            parameters = f"{config.y0};{config.x0};{config.A1};{config.t1};{config.A2};{config.t2};{config.A3};{config.t3};" \
                         f"{config.Min_gene_len};{config.Max_gene_len};{config.mean_gene_length};" \
                         f"{config.domen_length};{config.number_of_genes};{config.rotations_per_sample};" \
                         f"{config.minimal_expression};{config.maximal_expression}"
            file.write(parameters + '\n')
    except Exception as e:
        print(f"Error writing parameters to file: {e}")

def worker_func(args):
    """
    Функция для выполнения в пуле процессов.
    Аргументы:
    args: кортеж, содержащий индекс процесса j, путь к файлу file_path, и конфигурацию config
    """
    j, file_path, config = args  # Распаковываем переданные аргументы
    seed = int((time.time() * 1e6) % 1e6) + j
    np.random.seed(seed)
    random.seed(seed)
    if config.memorytype == 'buffer':
        chromosome = Modeler.Chromosome_generator(config.number_of_genes)
        z = one_chromosome_rework(chromosome = chromosome)
    elif config.memorytype == 'file':
        z = one_chromosome_rework(file_name = generated_data, chromnum = j)
    else:
        raise ValueError('unproper data in parametres file')
    append_number_to_file(z, file_path)
    print(f'Process {j} completed with z = {z}')
    return z


def run_in_processes_with_writer(num_tasks: int, file_path: str, config):
    """
    Функция для запуска процессов с помощью пула и записи результатов.
    Аргументы:
    num_tasks: количество задач (процессов)
    file_path: путь к файлу для записи
    config: объект конфигурации
    """
    args = [(j, file_path, config) for j in range(num_tasks)]
    if config.memorytype == 'buffer':
        pass
    elif config.memorytype == 'file':
        generated_data = Modeler.Generate_Chromosome_f(number_chroms = num_tasks, chrom_length = config.chromosome_length)
    with multiprocessing.Pool() as pool:
        results = pool.map(worker_func, args)
    print(f'All processes completed. Made {num_tasks} samples')
    return results


if __name__ == '__main__':
    config = Config()
    multiprocessing.freeze_support()
    timestamp = datetime.now().strftime('%y_%m_%d_%H_%M_%S') 
    file_path = f'out{timestamp}.txt'  # Убрано слово "output" из имени файла
    num_tasks = int(input('Start N: '))
    test_type = input('Choose test type (model/real, default is model): ') or 'model'  
    if test_type == 'model':
        pass  # Оставляем пустым, так как мы не используем arr_of_domens здесь
    elif test_type == 'real':
        # Предполагается, что вы реализуете метод load_real_data для загрузки реальных данных
        arr_of_domens = load_real_data()
    write_parameters_to_file(file_path, config)
    results = run_in_processes_with_writer(num_tasks, file_path, config)