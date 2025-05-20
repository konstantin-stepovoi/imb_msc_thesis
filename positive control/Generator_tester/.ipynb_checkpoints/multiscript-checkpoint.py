import os
import subprocess
from multiprocessing import Pool
from Generator import make_new_files_for_run
import time
from tqdm import tqdm
APP_PATH = "python3" 
from threading import Thread


def simulate_progress_for_generation(file_number):
    estimated_total_time = 6.0 * (file_number / 1000)  # 6 сек на 1000
    steps = 50  # количество "псевдошагов" прогресса
    sleep_time = estimated_total_time / steps

    for _ in tqdm(range(steps), desc="Generating data", bar_format="{l_bar}{bar}| {elapsed}", ncols=80):
        time.sleep(sleep_time)


def get_file_pairs(GENES_DIR, DOMAINS_DIR):
    genes_files = {file.split('.')[1] for file in os.listdir(GENES_DIR) if file.startswith("randomG.") and file.endswith(".bed")}
    domains_files = {file.split('.')[1] for file in os.listdir(DOMAINS_DIR) if file.startswith("randomD.") and file.endswith(".bed")}
    common_numbers = sorted(genes_files & domains_files)
    return common_numbers

def process_file_pair(file_number, iterations, GENES_DIR, DOMAINS_DIR):
    genes_file = os.path.join(GENES_DIR, f"randomG.{file_number}.bed")
    domains_file = os.path.join(DOMAINS_DIR, f"randomD.{file_number}.bed")

    try:
        process = subprocess.run(
            [APP_PATH, "onework.py"],
            input=f"{genes_file}\n{domains_file}\n{iterations}\n",
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if process.returncode == 0:
            print(f"N0 {file_number} Ok")
        else:
            print(f"Err on file {file_number}: {process.stderr}")
    except Exception as e:
        print(f"Err on file {file_number}: {e}")

def main():
    GENES_DIR = input('path to Genes directory: ')
    DOMAINS_DIR = input('path to Domains directory: ')
    file_number = int(input(f'number of realisations: '))
    MAP_TYPE = input('type of mapping (1-4, by default 2): ') or "2"
    MAP_TYPE = int(MAP_TYPE)
    percentage = float(int(input('Expression spread in % ')) / 100)

    print(f'\n Creating files for run...')

    gen_thread = Thread(target=make_new_files_for_run, args=(GENES_DIR, DOMAINS_DIR, file_number, percentage))
    gen_thread.start()
    simulate_progress_for_generation(file_number)
    gen_thread.join()
        
    iterations = int(input("Number of rotations: "))

    if not os.path.exists("result"):
        os.makedirs("result")

    file_numbers = get_file_pairs(GENES_DIR, DOMAINS_DIR)
    if not file_numbers:
        print("Err: no input files found, check the input dirs")
        return

    print(f"found {len(file_numbers)} file pairs to solve")

    num_cores = int(input('number of CPU cores available for this run: '))

    with open(f"metadata.txt", "w") as result_file:
        result_file.write("file_num\tmap_type\tpercentage\trotations\t\n")
        result_file.write(f'{file_number}\t{MAP_TYPE}\t{percentage}\t{iterations}')

    
    args = [(file_number, iterations, GENES_DIR, DOMAINS_DIR) for file_number in file_numbers]
    with Pool(processes=num_cores) as pool:
        pool.starmap(process_file_pair, args)

if __name__ == "__main__":
    main()
