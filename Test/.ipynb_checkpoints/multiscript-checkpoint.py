import os
import subprocess
import time
from multiprocessing import Pool, cpu_count

# Конфигурация
GENES_DIR = "Genes_Multiple"
DOMAINS_DIR = "Domains_Multiple"
APP_PATH = "python"  

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
    # Запрашиваем у пользователя тип маппинга
    print("Выберите тип маппинга:")
    print("1 - Ген относится к домену, в котором большая часть")
    print("2 - Экспрессия распределяется пропорционально длине в домене")
    print("3 - Экспрессия привязывается к точке старта")
    print("4 - Экспрессия привязывается к точке конца")
    
    while True:
        try:
            mapping_type = int(input("Введите номер метода (1-4): "))
            if mapping_type in {1, 2, 3, 4}:
                break
            else:
                print("Ошибка: Введите число от 1 до 4.")
        except ValueError:
            print("Ошибка: Введите корректное число.")
    
    # Записываем выбор пользователя в map_param.txt
    with open("map_param.txt", "w") as f:
        f.write(str(mapping_type))
    # Запрашиваем количество итераций у пользователя
    iterations = int(input("Введите количество итераций: "))
    
    if not os.path.exists("test"):
        os.makedirs("test")
    # Получаем список доступных пар файлов
    file_numbers = get_file_pairs()
    if not file_numbers:
        print("Не найдено пар файлов для обработки.")
        return

    print(f"Найдено {len(file_numbers)} пар файлов для обработки.")

    # Определяем количество доступных ядер
    num_cores = 6

    # Используем пул процессов для запуска одного экземпляра приложения на ядро
    with Pool(processes=num_cores) as pool:
        pool.starmap(process_file_pair, [(file_number, iterations) for file_number in file_numbers])

if __name__ == "__main__":
    start_time = time.time()
    main()
    print(f'ex time: {time.time() - start_time}')
