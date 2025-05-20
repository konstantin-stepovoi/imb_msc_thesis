# Interface.py
import numpy as np
import random
import time
import multiprocessing
from datetime import datetime
import Modeler
from parameters import Config, glob_config
from MathCore import one_chromosome_rework
from Parser import load_real_data  # Предположим, что этот модуль будет позже написан

config = Config()  # Глобальная конфигурация
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
    """
    j, file_path, config = args
    seed = int((time.time() * 1e6) % 1e6) + j
    np.random.seed(seed)
    random.seed(seed)

    z = one_chromosome_rework()  # Используем временную генерацию данных
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
        generated_data = Modeler.Generate_Chromosome_f(number_chroms=num_tasks, chrom_length=config.chromosome_length)
    with multiprocessing.Pool() as pool:
        results = pool.map(worker_func, args)
    print(f'All processes completed. Made {num_tasks} samples')
    return results

import os
output_folder = "lineartest"
os.makedirs(output_folder, exist_ok=True)


if __name__ == '__main__':
    config = Config()
    multiprocessing.freeze_support()
    num_tasks = 2000  # Количество задач на каждую итерацию
    rotations = [5000, 10000, 50000, 100000, 250000]

    for i in range(len(rotations)):
        glob_config.rotations_per_sample = rotations[i]  
        file_path = f'{output_folder}/Rotatest_{rotations[i]}.txt'  
        print(f"Starting test with rotor {rotations[i]}, file: {file_path}")
        
        # Записываем параметры в файл
        write_parameters_to_file(file_path, glob_config)
        
        # Запускаем процессы и записываем результаты в файл
        results = run_in_processes_with_writer(num_tasks, file_path, glob_config)
        
        # Инкремент глобального индекса для следующего теста
        print(f"Test with rotor {glob_config.rotations_per_sample} completed.")
