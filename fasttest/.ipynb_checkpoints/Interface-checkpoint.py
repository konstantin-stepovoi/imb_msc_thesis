import os
import multiprocessing
from JFL import many_chromosome_rework, process_bed_files  # Импортируем вашу функцию из JFL.py
import random
import numpy as np

def process_files(pair_index) -> None:
    # Генерация имен файлов на основе индекса пары
    gene_file = f"Genes/randomG.{pair_index:05d}.bed"
    domain_file = f"Domains/randomD.{pair_index:05d}.bed"
    try:
        # Получаем Domens и Genes с помощью функции из JFL.py
        Domens, Genes, expressions = process_bed_files(domain_file, gene_file)
        print('Splitted files')
        print(f'Domens: {len(Domens)}, {type(Domens), type(Domens[0])}')
        print(f'Genes: {len(Genes)}, {type(Genes), type(Genes[0])}')
        print(f'expresions: {len(expressions)}, {type(expressions), type(expressions[0])}')
        
        results = one_chromosome_rework(Domens, Genes, expressions, iterations)
        
        # Сохраняем результат в файл
        with open(f"results_many/{pair_index}.txt", "w") as result_file:
            result_file.write(f"z_crit, mu_orig, d_orig, mu_hor, d_hor, mu_vert, d_vert \n")
            for result in results:
                result_file.write(f"{result}\n")
        print(f'Processed file pair {pair_index}')

    except Exception as e:
        print(f"Error processing pair {pair_index}: {e}")

def main():
    # Папки, содержащие файлы с генами и доменами
    genes_dir = "Genes"
    domains_dir = "Domains"

    # Получаем все номера файлов (1-1000)
    files_count = len([name for name in os.listdir(genes_dir) if name.endswith('.bed')])

    # Создаем папку для результатов, если ее нет
    if not os.path.exists("results_differ"):
        os.makedirs("results_differ")

    # Используем multiprocessing для параллельной обработки
    # Получаем метод запуска multiprocessing
    method = multiprocessing.get_start_method()

    # Обработка всех пар файлов в параллели
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(process_files, range(1, files_count + 1))  # Обрабатываем все файлы

if __name__ == "__main__":
    multiprocessing.freeze_support()  # Поддержка для заморозки в Windows
    main()
