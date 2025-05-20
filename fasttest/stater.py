import os

# Папка с файлами
directory = '.'

# Массив для хранения первого числа из третьей строки
numbers = []

# Пробегаемся по всем .txt файлам в папке
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        filepath = os.path.join(directory, filename)
        with open(filepath, 'r') as file:
            lines = file.readlines()
            if len(lines) >= 4:
                # Третья строка
                line = lines[3]
                # Убираем лишние символы и преобразуем в список чисел
                line = line.strip().lstrip('[').rstrip(']').split(',')
                try:
                    # Добавляем первое число в массив
                    first_number = float(line[0])
                    numbers.append(first_number)
                except (ValueError, IndexError):
                    print(f"Ошибка обработки строки в файле {filename}: {line}")

# Сохраняем числа в stats.txt
with open('stats.txt', 'w') as stats_file:
    for number in numbers:
        stats_file.write(f"{number}\n")
