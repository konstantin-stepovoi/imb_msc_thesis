#!/usr/bin/python3
import os
from Gene_rotor_matematics import Weilford_function, process_bed_files  
import random
import time
import numpy as np

def process_files(name1, name2, iterations, result_dir) -> None:
    q = int(name1.split('.')[1])
    Domens, Genes, expressions = process_bed_files(name2, name1)
    results = Weilford_function(Domens, Genes, expressions, iterations)

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
    result_dir = input('result_dir: ').strip()   # <- отличие от прошлого!
    process_files(name1, name2, iterations, result_dir)

if __name__ == "__main__":
    main()
