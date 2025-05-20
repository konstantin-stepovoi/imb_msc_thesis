import os
import numpy as np
import matplotlib.pyplot as plt

def process_txt_files():
    # Собираем float числа из всех .txt файлов
    numbers = []

    for filename in os.listdir('.'):  # Перебираем файлы в текущей директории
        if filename.endswith('.txt'):
            try:
                with open(filename, 'r') as file:
                    number = float(file.read().strip())  # Читаем и конвертируем в float
                    numbers.append(number)
            except ValueError:
                print(f"Файл {filename} не содержит валидного float числа и будет пропущен.")

    if not numbers:
        print("Нет валидных данных для обработки.")
        return

    # Вычисляем среднее и дисперсию
    mean = np.mean(numbers)
    variance = np.var(numbers)
    print(f"Среднее: {mean}")
    print(f"Дисперсия: {variance}")

    # Строим гистограмму
    plt.figure(figsize=(10, 6))
    plt.hist(numbers, bins=30, color='blue', edgecolor='black', alpha=0.7)
    plt.title("Гистограмма значений", fontsize=16)
    plt.xlabel("Значение", fontsize=14)
    plt.ylabel("Частота", fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Сохраняем гистограмму как PNG
    plt.savefig("histogram.png")
    print("Гистограмма сохранена в файл 'histogram.png'.")

if __name__ == "__main__":
    process_txt_files()