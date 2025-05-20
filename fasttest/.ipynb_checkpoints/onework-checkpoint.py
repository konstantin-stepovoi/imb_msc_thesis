import os
from JFL2 import one_chromosome_rework_new, process_bed_files  # Импортируем вашу функцию из JFL.py
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
        print(f"Results: {results}")  # Проверяем, что вернула функция

      
        ###############################
        # Сохраняем результат в файл
        with open(f"prefiltered_1000/{q:05d}.txt", "w") as result_file:
            result_file.write(f'iterations: {iter} \n')
            result_file.write(f"z_crit, mu_orig, d_orig, mu_hor, d_hor, mu_vert, d_vert \n")
            for result in results:
                result_file.write(f"{result}\n")
        print(f'Processed file pair {q}')
    except Exception as e:
        print(f"Error processing {name1}: {e}")

def main():
    if not os.path.exists("prefiltered_1000"):
        os.makedirs("prefiltered_1000")
        
    name1 = input(f'Genes_file_path:')
    name2 = input(f'Domains_file_path:')
    iter = int(input(f'number_iterations'))
    start = time.time()
    process_files(name1, name2, iter)
    print(time.time()-start)
if __name__ == "__main__":
    main()

