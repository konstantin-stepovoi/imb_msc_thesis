import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import anderson, norm

# Получаем список всех txt файлов в текущей папке
filenames = [f for f in os.listdir('.') if f.endswith('.txt')]

# Массив для хранения значений
values = []
# Проходим по файлам и извлекаем нужные данные
for filename in filenames:
    with open(filename, 'r') as file:
        lines = file.readlines()
        if len(lines) > 1:
            try:
                value = float(lines[1].split()[0])  # Берем первое число из второй строки
                values.append(value)
            except ValueError:
                print(f"Ошибка чтения числа в файле {filename}")

# Преобразуем в numpy массив
values.append(1.95)
values = np.array(values)
print(values)
# Вычисляем статистику
mean_value = np.mean(values)
std_value = np.std(values, ddof=1)  # Выборочное стандартное отклонение
variance = np.var(values, ddof=1)

# Проверка на нормальность методом Андерсона-Дарлинга
result = anderson(values)

# Построение гистограммы
plt.figure(figsize=(10, 6))
bins = 45
count, bins_edges, _ = plt.hist(values, bins=bins, color='skyblue', edgecolor='black', density=True, label='Histogram')

# Фит нормального распределения
x = np.linspace(min(values), max(values), 1000)
gaussian = norm.pdf(x, mean_value, std_value)
plt.plot(x, gaussian, 'k--', linewidth=1.5, label='Normal Fit')

# Оформление
percentage = int(input(f'whats the percentage: '))
plt.title(f"Разброс при генерации: {percentage}%\n"
          f"Mean: {mean_value:.2f}, Variance: {variance:.4f}, Std: {std_value:.4f}\n"
          f"Anderson-Darling Statistic: {result.statistic:.2f}", fontsize=14)
plt.xlabel("Значение Z (Unpadded)")
plt.ylabel("Плотность вероятности")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# Сохраняем график в PNG с высоким разрешением
plt.savefig("histogram.png", dpi=300, bbox_inches='tight')

# Показываем график
plt.show()

# Вывод результатов теста Андерсона-Дарлинга
print("Anderson-Darling:")
for i in range(len(result.critical_values)):
    significance_level = result.significance_level[i]
    critical_value = result.critical_values[i]
    #print(f"Significance Level {significance_level}%: Critical Value = {critical_value:.3f}")
if result.statistic < result.critical_values[2]:  # 5% уровень значимости
    #print("Распределение можно считать нормальным на уровне значимости 5%.")
else:
    print("Распределение значительно отклоняется от нормального на уровне значимости 5%.")

# Подсчёт статистик
n_total  = len(values)
n_neg = np.sum(values < 0)
n_less_1_96 = np.sum(values < 1.96)
n_greater_2_58 = np.sum(values > 2.58)
print(n_total)
# Вывод результатов
print(f"Доля значений < 1.96: {n_less_1_96} ({n_less_1_96 / n_total:.4f})")
print(f"Доля значений > 2.58: {n_greater_2_58} ({n_greater_2_58 / n_total:.4f})")
print(f"Общая доля: {n_greater_2_58 + n_less_1_96} ({(n_greater_2_58 + n_less_1_96) / n_total:.4f})")
