#!/usr/bin/python3
import os
from Gene_rotor_matematics import run_with_config, process_bed_files
import random
import time
import numpy as np

def process_files(name1, name2, iterations, result_dir, mapper_method) -> None:
    q = int(name1.split('.')[1])
    Domens, Genes, expressions = process_bed_files(name2, name1)
    
    # Запускаем через обёртку с выбранным методом
    results = run_with_config(Domens, Genes, expressions, iterations, mapper_method=mapper_method)

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    result_file_path = os.path.join(result_dir, f"{q:05d}.txt")
    with open(result_file_path, "w") as result_file:
        result_file.write("z_crit\tmu_orig\td_orig\tmu_hor\td_hor\tresult\tcheckstat\n")
        result_file.write("\t".join(map(str, results)) + "\n")
    print(f'Processed file pair {q}')


def main():
    name1 = input('Genes_file_path: ').strip()
    name2 = input('Domains_file_path: ').strip()
    iterations = int(input('number_iterations: ').strip())
    result_dir = input('result_dir: ').strip()
    mapper_method = input('mapper_method: ').strip()
    process_files(name1, name2, iterations, result_dir, mapper_method)

def temp():
    name1 = input('Genes_file_path: ').strip()
    name2 = input('Domains_file_path: ').strip()
    iterations = int(input('number_iterations: ').strip())
    result_dir = input('result_dir: ').strip()
    mapper_method = input('mapper_method: ').strip()
    process_files(name1, name2, iterations, result_dir, mapper_method)
    
if __name__ == "__main__":
    main()
