#!/usr/bin/python3
# Импорты
# Interface.py
import numpy as np
import random, time, multiprocessing
from multiprocessing.managers import ListProxy
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

def worker_func(j: int, result_list: ListProxy, file_path: str, config: Config):
    seed = int((time.time() * 1e6) % 1e6) + j
    np.random.seed(seed)
    random.seed(seed)
    
    # Генерация новой хромосомы для каждого процесса
    chromosome = Modeler.Chromosome_generator(config.number_of_genes)
    
    # Теперь передаем chromosome в функцию
    z = one_chromosome_rework(chromosome)
    append_number_to_file(z, file_path)
    result_list.append(z)
    print(f'Process {j} completed with z = {z}')

def run_in_processes_with_writer(num_tasks: int, file_path: str, config: Config):
    manager = multiprocessing.Manager()
    results = manager.list()
    
    processes = []
    for j in range(num_tasks):
        p = multiprocessing.Process(target=worker_func, args=(j, results, file_path, config))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

    print(f'All processes completed. made {num_tasks} samples')

if __name__ == '__main__':
    config = Config()
    multiprocessing.freeze_support()

    # Генерация имени файла на основе текущей даты и времени
    timestamp = datetime.now().strftime('%y_%m_%d_%H_%M_%S')  # Изменен формат временной метки
    file_path = f'{timestamp}.txt'  # Убрано слово "output" из имени файла

    num_tasks = int(input('Start N: '))
    test_type = input('Choose test type (model/real, default is model): ') or 'model'  # Установлено значение по умолчанию

    if test_type == 'model':
        # Здесь вызов функции Modeler для генерации данных будет происходить внутри worker_func
        pass  # Можно оставить пустым, так как мы не используем arr_of_domens здесь
    elif test_type == 'real':
        arr_of_domens = load_real_data()  # Предполагается, что вы реализуете этот метод в Parser

    write_parameters_to_file(file_path, config)  # Запись параметров в файл
    run_in_processes_with_writer(num_tasks, file_path, config)  # Передаем config вместо arr_of_domens
