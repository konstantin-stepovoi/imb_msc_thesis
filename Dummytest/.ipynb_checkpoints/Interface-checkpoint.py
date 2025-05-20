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

    # Устанавливаем seed для каждого процесса на основе текущего времени и индекса j
    seed = int((time.time() * 1e6) % 1e6) + j
    np.random.seed(seed)
    random.seed(seed)

    # Генерация новой хромосомы для каждого процесса
    chromosome = Modeler.Chromosome_generator(config.number_of_genes)

    # Теперь передаем chromosome в функцию для обработки
    z = one_chromosome_rework(chromosome)

    # Запись результата в файл
    append_number_to_file(z, file_path)

    # Печатаем, что процесс завершился
    print(f'Process {j} completed with z = {z}')
    
    # Возвращаем результат для последующего использования
    return z


def run_in_processes_with_writer(num_tasks: int, file_path: str, config):
    """
    Функция для запуска процессов с помощью пула и записи результатов.
    Аргументы:
    num_tasks: количество задач (процессов)
    file_path: путь к файлу для записи
    config: объект конфигурации
    """
    # Подготавливаем аргументы для каждого процесса
    args = [(j, file_path, config) for j in range(num_tasks)]

    # Используем Pool для запуска рабочих процессов
    with multiprocessing.Pool() as pool:
        # Передаем задачи в пул процессов и собираем результаты
        results = pool.map(worker_func, args)

    # Когда все процессы завершены, результаты уже записаны
    print(f'All processes completed. Made {num_tasks} samples')

    return results


if __name__ == '__main__':
    # Загружаем конфигурацию
    config = Config()

    # Необходимая настройка для работы с замороженными объектами в Windows
    multiprocessing.freeze_support()

    # Генерация имени файла на основе текущей даты и времени в формате гг_мм_дд_чч_мм_сс
    timestamp = datetime.now().strftime('%y_%m_%d_%H_%M_%S')  # Изменен формат временной метки
    file_path = f'out{timestamp}.txt'  # Убрано слово "output" из имени файла

    # Получаем количество задач от пользователя
    num_tasks = int(input('Start N: '))

    # Получаем тип теста от пользователя (model/real)
    test_type = input('Choose test type (model/real, default is model): ') or 'model'  # Установлено значение по умолчанию

    # Определяем действия в зависимости от типа теста
    if test_type == 'model':
        # Здесь вызов функции Modeler для генерации данных будет происходить внутри worker_func
        pass  # Оставляем пустым, так как мы не используем arr_of_domens здесь
    elif test_type == 'real':
        # Предполагается, что вы реализуете метод load_real_data для загрузки реальных данных
        arr_of_domens = load_real_data()

    # Записываем параметры в файл
    write_parameters_to_file(file_path, config)

    # Запуск процессов с записью результатов
    results = run_in_processes_with_writer(num_tasks, file_path, config)