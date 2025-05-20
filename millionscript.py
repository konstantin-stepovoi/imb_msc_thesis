import os
import subprocess
from multiprocessing import Pool

# Конфигурация
GENES_DIR = "Genes"
DOMAINS_DIR = "Domains"
APP_PATH = "python3"  # Путь к вашему приложению, если это python-скрипт
OUTPUT_DIR = "milliontest_1K"


def get_file_pairs():
    """
    Сканирует папки GENES_DIR и DOMAINS_DIR и формирует список доступных пар файлов.
    Возвращает список кортежей с номерами файлов.
    """
    genes_files = {file.split('.')[1] for file in os.listdir(GENES_DIR) if file.startswith("randomG.") and file.endswith(".bed")}
    domains_files = {file.split('.')[1] for file in os.listdir(DOMAINS_DIR) if file.startswith("randomD.") and file.endswith(".bed")}
    common_numbers = sorted(genes_files & domains_files)
    return common_numbers


def process_file_pair(file_number, iterations):
    """
    Запускает консольное приложение для обработки одной пары файлов 1000 раз.
    """
    genes_file = os.path.join(GENES_DIR, f"randomG.{file_number}.bed")
    domains_file = os.path.join(DOMAINS_DIR, f"randomD.{file_number}.bed")
    output_file = os.path.join(OUTPUT_DIR, f"output_{file_number}.txt")

    try:
        #with open(output_file, "a") as out_file:
            #for i in range(1): #Пока это вручную, потом проверю и добавлю 1000
                
        process = subprocess.run(
                [APP_PATH, "onework.py", genes_file, domains_file, str(50)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        if process.returncode == 0:
            print(f'Done file pair {file_number}')
                    #out_file.write(process.stdout.strip() + "\n")
        else:
            print(f"Ошибка при обработке {file_number}, итерация {i + 1}: {process.stderr}")

    

    except Exception as e:
        print(f"Ошибка при запуске обработки для {file_number}: {e}")


def main():
    # Запрашиваем количество итераций у пользователя
    iterations = 20 
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Получаем список доступных пар файлов
    file_numbers = get_file_pairs()
    if not file_numbers:
        print("Не найдено пар файлов для обработки.")
        return

    print(f"Найдено {len(file_numbers)} пар файлов для обработки.")

    # Определяем количество доступных ядер
    num_cores = int(input(f'What is current number of cores: '))

    # Используем пул процессов для параллельного запуска
    with Pool(processes=num_cores) as pool:
        pool.starmap(process_file_pair, [(file_number, iterations) for file_number in file_numbers])


if __name__ == "__main__":
    main()
