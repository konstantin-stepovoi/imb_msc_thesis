#!/usr/bin/python.exe
#!/usr/bin/env python3
# Проверка необходимых библиотек и python executable


required_modules = ['os', 'sys', 'subprocess', 'time', 'json', 'threading',
                    'multiprocessing', 'datetime', 'tqdm', 'glob', 'random', 'numpy']

missing = []
for module in required_modules:
    try:
        __import__(module)
    except ImportError:
        missing.append(module)

if missing:
    print(f"Moduled reqired for run: {', '.join(missing)} are not installed")
    print("Try to install them with:")
    print(f"   pip install {' '.join(missing)}")
    sys.exit(1)

import sys
import os
import subprocess
import time
import json
import threading
from multiprocessing import Pool
from datetime import datetime
from tqdm import tqdm
import glob
import random
import numpy as np
import logging
from JFL3 import one_chromosome_rework_new, process_bed_files
import shutil

print(f"Python executable: {sys.executable}")

# Настройка логирования
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[logging.StreamHandler(), logging.FileHandler('process_log.txt', mode='w')])


SETTINGS_FILE = "settings.json"
DEFAULT_SETTINGS = {
    "genes_output_dir": "",
    "domains_output_dir": "",
    "genes_input_dir": "",
    "domains_input_dir": "",
    "mapper_method": "proportional",
    "results_dir": "result",
    "analysis_name": "default_analysis",
    "cpu_cores": 2,
    "percentage": 0.05,
    "iterations": 10,
    "analysis_type": "genome-wide",  # Добавили новую настройку
}

APP_PATH = sys.executable

def load_settings():
    if os.path.exists(SETTINGS_FILE):
        with open(SETTINGS_FILE, "r") as f:
            return json.load(f)
    return DEFAULT_SETTINGS.copy()

def save_settings(settings):
    with open(SETTINGS_FILE, "w") as f:
        json.dump(settings, f, indent=4)

def print_help():
    print("Available commands:")
    print("  show-settings               Show current settings")
    print("  set-[param] [value]        Set a configuration parameter")
    print("     Available params:")
    print("       genes_output_dir, domains_output_dir")
    print("       genes_input_dir, domains_input_dir")
    print("       mapper_method [middle|proporcional|start|stop]")
    print("       results_dir, analysis_name")
    print("       cpu_cores, percentage, iterations")
    print("       analysis_type [genome-wide|chromosome-wide]")  # Добавили вывод для analysis_type
    print("  generate                    Start data generation")
    print("  analyze                     Start analysis")
    print("  exit                        Exit the program")

def is_absolute(path):
    return os.path.isabs(path)

def simulate_progress_for_generation(file_number):
    estimated_total_time = 6.0 * (file_number / 1000)
    steps = 50
    sleep_time = estimated_total_time / steps
    for _ in tqdm(range(steps), desc="Generating data", bar_format="{l_bar}{bar}| {elapsed}", ncols=80):
        time.sleep(sleep_time)

def get_file_pairs(genes_dir, domains_dir):
    genes_files = {file.split('.')[1] for file in os.listdir(genes_dir) if file.startswith("randomG.") and file.endswith(".bed")}
    domains_files = {file.split('.')[1] for file in os.listdir(domains_dir) if file.startswith("randomD.") and file.endswith(".bed")}
    return sorted(genes_files & domains_files)

def process_file_pair(file_number, iterations, genes_dir, domains_dir, results_dir):
    genes_file = os.path.join(genes_dir, f"randomG.{file_number}.bed")
    domains_file = os.path.join(domains_dir, f"randomD.{file_number}.bed")
    
    # Проверка существования файлов
    if not os.path.exists(genes_file):
        logging.error(f"Gene file {genes_file} does not exist")
        return
    if not os.path.exists(domains_file):
        logging.error(f"Domain file {domains_file} does not exist")
        return
    
    # Логирование пути к файлам и передаваемых данных
    logging.info(f"Processing file pair {file_number}")
    logging.info(f"Genes file: {genes_file}")
    logging.info(f"Domains file: {domains_file}")
    logging.info(f"Iterations: {iterations}")
    logging.info(f"Results directory: {results_dir}")
    
    try:
        process = subprocess.run(
            [APP_PATH, "onework.py"],
            input=f"{genes_file}\n{domains_file}\n{iterations}\n{results_dir}\n",
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        
        # Логирование статуса выполнения
        if process.returncode == 0:
            logging.info(f"File pair {file_number} processed successfully")
        else:
            logging.error(f"Error processing file pair {file_number}: {process.stderr}")
    except Exception as e:
        logging.error(f"Exception occurred while processing file pair {file_number}: {e}")

def run_generator(settings):
    from Generator import make_new_files_for_run
    file_number = int(input("Number of realisations to generate: "))
    settings['percentage'] = float(input("Expression spread in %: ")) / 100
    save_settings(settings)

    gen_thread = threading.Thread(target=make_new_files_for_run, args=(
        settings['genes_output_dir'],
        settings['domains_output_dir'],
        file_number,
        settings['percentage']))
    gen_thread.start()
    simulate_progress_for_generation(file_number)
    gen_thread.join()

def run_analysis(settings):
    file_numbers = get_file_pairs(settings['genes_input_dir'], settings['domains_input_dir'])
    if not file_numbers:
        print("Err: no input files found, check the input dirs")
        return
        
    print(f'found {len(file_numbers)} to analyze')

    if not os.path.exists(settings['results_dir']):
        os.makedirs(settings['results_dir'])

    with open("cur_metadata.txt", "w") as meta:
        meta.write(f"Analysis: {settings['analysis_name']}\n")
        meta.write(f"Datetime: {datetime.now()}\n")
        for key, val in settings.items():
            meta.write(f"{key}: {val}\n")

    args = [(fnum, settings['iterations'], settings['genes_input_dir'], settings['domains_input_dir'], settings['results_dir']) for fnum in file_numbers]
    with Pool(processes=settings['cpu_cores']) as pool:
        pool.starmap(process_file_pair, args)

def collect_bcp_data(results_dir):
    bcp_output_dir = os.path.join(os.getcwd(), "bcp_dt")
    os.makedirs(bcp_output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(bcp_output_dir, f"{timestamp}.txt")

    result_files = sorted(glob.glob(os.path.join(results_dir, "*.txt")))

    with open(output_file, "w") as out:
        for file in result_files:
            try:
                with open(file, "r") as f:
                    lines = f.readlines()
                    if len(lines) > 1:
                        out.write(lines[1].strip() + "\n")
            except Exception as e:
                print(f"Could not read from {file}: {e}")

    print(f"BCP summary written to: {output_file}")


def chromosome_wide_analysis(settings):
    genes_dir = os.path.abspath(settings['input_genes_dir'])
    domains_dir = os.path.abspath(settings['input_domains_dir'])
    results_dir = os.path.abspath(settings['results_dir'])
    temp_genes_dir = os.path.abspath('Genes_temp')
    temp_domains_dir = os.path.abspath('Domains_temp')

    # создаём временные директории
    os.makedirs(temp_genes_dir, exist_ok=True)
    os.makedirs(temp_domains_dir, exist_ok=True)

    # находим уникальные хромосомы по первому gene-файлу
    sample_file = next((f for f in os.listdir(genes_dir) if os.path.isfile(os.path.join(genes_dir, f))), None)
    if not sample_file:
        raise RuntimeError("No gene files found in input_genes_dir")

    chrom_set = set()
    with open(os.path.join(genes_dir, sample_file)) as f:
        for line in f:
            if line.strip():
                chrom = line.split()[0]
                chrom_set.add(chrom)

    for chrom in sorted(chrom_set):
        print(f"--- Processing {chrom} ---")

        # создаём results/chrN
        chrom_result_dir = os.path.join(results_dir, chrom)
        os.makedirs(chrom_result_dir, exist_ok=True)

        # чистим временные папки
        for folder in [temp_genes_dir, temp_domains_dir]:
            for file in os.listdir(folder):
                os.remove(os.path.join(folder, file))

        # фильтруем и копируем только нужные строки по хромосоме
        for file in os.listdir(genes_dir):
            in_path = os.path.join(genes_dir, file)
            out_path = os.path.join(temp_genes_dir, file)
            with open(in_path) as src, open(out_path, 'w') as dst:
                for line in src:
                    if line.startswith(chrom + '\t'):
                        dst.write(line)

        for file in os.listdir(domains_dir):
            in_path = os.path.join(domains_dir, file)
            out_path = os.path.join(temp_domains_dir, file)
            with open(in_path) as src, open(out_path, 'w') as dst:
                for line in src:
                    if line.startswith(chrom + '\t'):
                        dst.write(line)

        # запускаем анализ на одну хромосому
        chrom_settings = settings.copy()
        chrom_settings['input_genes_dir'] = temp_genes_dir
        chrom_settings['input_domains_dir'] = temp_domains_dir
        chrom_settings['results_dir'] = chrom_result_dir

        run_analysis(chrom_settings)

    # удаляем временные директории
    shutil.rmtree(temp_genes_dir)
    shutil.rmtree(temp_domains_dir)

    print("Chromosome-wide analysis completed.")


def main():
    print("*******************")
    print("GeneRotor is ready.")
    print("*******************")
    print("Type \"--help\" to see options")

    settings = load_settings()

    while True:
        cmd = input(">> ").strip()
        if cmd == "--help":
            print_help()
        elif cmd == "show-settings":
            for k, v in settings.items():
                print(f"{k}: {v}")
        elif cmd.startswith("set-"):
            try:
                param, value = cmd[4:].split(None, 1)
                if 'dir' in param and not is_absolute(value):
                    print("Please provide an absolute path for directories.")
                    continue
                if param == 'cpu_cores' or param == 'iterations':
                    settings[param] = int(value)
                elif param == 'percentage':
                    settings[param] = float(value)
                elif param == 'mapper_method' and value in ['middle', 'proporcional', 'start', 'stop']:
                    settings[param] = value
                elif param == 'analysis_type' and value in ['genome-wide', 'chromosome-wide']:  # Проверка значения для analysis_type
                    settings[param] = value
                elif param in settings:
                    settings[param] = value
                else:
                    print("Unknown parameter")
                    continue
                save_settings(settings)
            except ValueError:
                print("Usage: set-[param] [value]")
        elif cmd == "generate":
            run_generator(settings)
        elif cmd == "analyze":
            if settings['analysis_type'] == 'genome-wide':
                run_analysis(settings)
                collect_bcp_data(settings['results_dir'])
            elif settings['analysis_type'] == 'chromosome-wide':
                chromosome_wide_analysis(settings)
            else:
                raise ValueError("Unknown analysis_type in settings")
        elif cmd == "exit":
            break
        else:
            print("Unknown command. Type --help for options.")

if __name__ == "__main__":
    main()
