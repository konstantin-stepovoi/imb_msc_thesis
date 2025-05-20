#!/usr/bin/python3
import os
from JFL3 import one_chromosome_rework_new, process_bed_files  
import random
import time
import numpy as np

def process_files(name1, name2, iter) -> None:
    gene_file = name1
    domain_file = name2
    iterations = iter
    q = int(gene_file.split('.')[1])
    try:
        # Получаем Domens и Genes с помощью функции из JFL.py
        Domens, Genes, expressions = process_bed_files(domain_file, gene_file)
        print('Splitted files')
        print(f'Domens: {len(Domens)}, {type(Domens), type(Domens[0])}')
        print(f'Genes: {len(Genes)}, {type(Genes), type(Genes[0])}')
        print(f'expresions: {len(expressions)}, {type(expressions), type(expressions[0])}')
        results = one_chromosome_rework_new(Domens, Genes, expressions, iterations)

        # Сохраняем результат в файл
        with open(f"result/{q:05d}.txt", "w") as result_file:
            result_file.write("z_crit\tmu_orig\td_orig\tmu_hor\td_hor\tmu_vert\td_vert\n")
            result_file.write("\t".join(map(str, results)) + "\n")
        print(f'Processed file pair {q}')
    except Exception as e:
        print(f"Error processing {name1}: {e}")

def main():
    if not os.path.exists("result"):
        os.makedirs("result")
        
    name1 = input(f'Genes_file_path:')
    name2 = input(f'Domains_file_path:')
    iter = int(input(f'number_iterations'))
    process_files(name1, name2, iter)
if __name__ == "__main__":
    main()

