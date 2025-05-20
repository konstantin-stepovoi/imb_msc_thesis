import os
import subprocess
from multiprocessing import Pool, cpu_count

# Конфигурация
GENES_DIR = "Genes"
DOMAINS_DIR = "Domains"
APP_PATH = "python"  # Путь к вашему приложению, если это python-скрипт

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
    Запускает консольное приложение для обработки одной пары файлов.
    """
    genes_file = os.path.join(GENES_DIR, f"randomG.{file_number}.bed")
    domains_file = os.path.join(DOMAINS_DIR, f"randomD.{file_number}.bed")

    try:
        # Запускаем консольное приложение и передаем аргументы
        process = subprocess.run(
            [APP_PATH, "onework.py"],  # Замените "onework.py" на имя вашего скрипта
            input=f"{genes_file}\n{domains_file}\n{iterations}\n",  # Передача данных в stdin
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if process.returncode == 0:
            print(f"Файлы {file_number} обработаны успешно.")
        else:
            print(f"Ошибка при обработке файлов {file_number}: {process.stderr}")
    except Exception as e:
        print(f"Ошибка при запуске обработки для {file_number}: {e}")

def main():
    # Запрашиваем количество итераций у пользователя
    iterations = int(input("Введите количество итераций: "))

    # Получаем список доступных пар файлов
    file_numbers = get_file_pairs()
    if not file_numbers:
        print("Не найдено пар файлов для обработки.")
        return

    print(f"Найдено {len(file_numbers)} пар файлов для обработки.")

    # Определяем количество доступных ядер
    num_cores = cpu_count() - 2

    # Используем пул процессов для запуска одного экземпляра приложения на ядро
    with Pool(processes=num_cores) as pool:
        pool.starmap(process_file_pair, [(file_number, iterations) for file_number in file_numbers])

if __name__ == "__main__":
    main()
