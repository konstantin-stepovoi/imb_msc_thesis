import random
import numpy as np
import multiprocessing
import time

def process_bed_files(file1, file2):
    # Открываем файлы и читаем строки
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    # Создаем два пустых списка
    list1 = []
    list2 = []
    list3 = []

    # Функция для обработки строк и заполнения списка
    def process_lines(lines, result_list):
        for line in lines:
            if line.startswith("chr"):
                parts = line.split()
                if len(parts) >= 3:
                    # Получаем chr3 и два числа
                    chr_part = parts[0]
                    try:
                        num1 = int(parts[1])
                        num2 = int(parts[2])
                        result_list.append((num1, num2))
                    except ValueError:
                        continue

    def process_lines2(lines, result_list):
        for line in lines:
            if line.startswith("chr"):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        num3 = float(parts[3])
                        result_list.append(num3)
                    except ValueError:
                        continue
    # Обрабатываем оба файла
    process_lines(lines1, list1)
    process_lines(lines2, list2)
    process_lines2(lines2, list3)
    return list1, list2, list3

def map_genes_to_domains2(Domens, Genes, expressions):
    # Преобразуем входные данные в NumPy массивы
    Domens = np.array(Domens)
    Genes = np.array(Genes)
    expressions = np.array(expressions)

    # Создаем массив для хранения доменов
    arr_of_domens = np.empty(len(Domens), dtype=object)

    # Проходим по каждому домену
    for domain_index, (dom_start, dom_end) in enumerate(Domens):
        # Находим индексы генов, которые пересекаются с доменом
        overlaps = (Genes[:, 1] >= dom_start) & (Genes[:, 0] <= dom_end)

        # Извлекаем данные о пересекающихся генах
        overlapping_genes = Genes[overlaps]
        overlapping_expressions = expressions[overlaps]

        this_domen = []

        # Обрабатываем пересекающиеся гены
        for (gene_start, gene_end), expression in zip(overlapping_genes, overlapping_expressions):
            gene_length = gene_end - gene_start

            # Полное вхождение гена в домен
            if dom_start <= gene_start and gene_end <= dom_end:
                this_domen.append(expression)

            # Полное покрытие домена геном
            elif gene_start <= dom_start and dom_end <= gene_end:
                overlap_length = dom_end - dom_start
                this_domen.append(expression * (overlap_length / gene_length))

            # Ген перекрывает левую границу домена
            elif gene_start < dom_start < gene_end <= dom_end:
                overlap_length = gene_end - dom_start
                this_domen.append(expression * (overlap_length / gene_length))

            # Ген перекрывает правую границу домена
            elif dom_start <= gene_start < dom_end < gene_end:
                overlap_length = dom_end - gene_start
                this_domen.append(expression * (overlap_length / gene_length))

        arr_of_domens[domain_index] = np.array(this_domen)

    return arr_of_domens

def map_genes_to_domains(Domens, Genes, expressions):
    # Создаем пустой массив для доменов
    arr_of_domens = np.empty(len(Domens), dtype=object)

    # Указатель на текущий ген
    gene_index = 0

    # Проходим по каждому домену
    for domain_index, (dom_start, dom_end) in enumerate(Domens):
        this_domen = []

        # Сдвигаем указатель генов, пока ген полностью "до" домена
        while gene_index < len(Genes) and Genes[gene_index][1] < dom_start:
            gene_index += 1

        # Обрабатываем гены, которые пересекаются с доменом
        temp_index = gene_index  # Сохраняем текущую позицию указателя
        while temp_index < len(Genes) and Genes[temp_index][0] < dom_end:
            gene_start, gene_end = Genes[temp_index]
            gene_length = gene_end - gene_start
            expression = expressions[temp_index]

            # Полное вхождение гена в домен
            if dom_start <= gene_start and gene_end <= dom_end:
                this_domen.append(expression)

            # Полное покрытие домена геном
            elif gene_start <= dom_start and dom_end <= gene_end:
                overlap_length = dom_end - dom_start
                this_domen.append(expression * (overlap_length / gene_length))

            # Ген перекрывает левую границу домена
            elif gene_start < dom_start < gene_end <= dom_end:
                overlap_length = gene_end - dom_start
                this_domen.append(expression * (overlap_length / gene_length))

            # Ген перекрывает правую границу домена
            elif dom_start <= gene_start < dom_end < gene_end:
                overlap_length = dom_end - gene_start
                this_domen.append(expression * (overlap_length / gene_length))

            # Переход к следующему гену
            temp_index += 1

        arr_of_domens[domain_index] = np.array(this_domen)

    return arr_of_domens

def count_statvalues_horis2(arr_of_domens: np.ndarray) -> tuple:
    # Фильтруем массивы, убирая None и пустые массивы
    valid_domens = [np.array(domen) for domen in arr_of_domens if domen is not None and len(domen) > 0]
    
    if not valid_domens:  # Проверяем, остались ли валидные массивы
        return None, None
    
    # Считаем средние и дисперсии для каждого домена с помощью list comprehension
    means = [np.mean(domen) for domen in valid_domens]
    variances = [np.var(domen) for domen in valid_domens]
    
    # Вычисляем итоговые значения
    overall_mean = np.mean(means)
    overall_dispersion = np.mean(variances)

    return overall_mean, overall_dispersion

def pad_and_stack(arr_of_domens):
    max_length = max(len(domen) for domen in arr_of_domens)
    padded_array = np.array([np.pad(domen, (0, max_length - len(domen)), constant_values=np.nan)
                             for domen in arr_of_domens])
    return padded_array
    
def count_statvalues_horis(arr_of_domens: np.ndarray) -> tuple:
    means = np.nanmean(arr_of_domens, axis=1)
    variances = np.nanvar(arr_of_domens, axis=1)

    overall_mean = np.nanmean(means)
    overall_dispersion = np.nanmean(variances)
    return overall_mean, overall_dispersion

def z_criterium(mu1: float, mu2: float, sig1: float, sig2: float, n: int) -> float:
    return (mu1 - mu2) / np.sqrt(sig1/n + sig2/n)

def rotate_genes(Domens, Genes, expressions, shift_value):
    # Находим общую длину хромосомы (конец последнего гена)
    chrom_length = max(gene[1] for gene in Genes)

    # Сдвигаем координаты всех генов
    new_genes = []
    for start, end in Genes:
        new_start = (start + shift_value) % chrom_length
        new_end = (end + shift_value) % chrom_length
        new_genes.append((new_start, new_end))

    # Массив expressions должен быть также сдвинут
    new_expressions = expressions[-shift_value % len(expressions):] + expressions[:-shift_value % len(expressions)]

    # Вызываем функцию маппинга с новыми сдвинутыми данными
    arr_of_domens = map_genes_to_domains(Domens, new_genes, new_expressions)
    filtered_domens = pad_and_stack(arr_of_domens)
    return filtered_domens

#ОБИДНО: если бы домены были одинаковой длины, я могбы это
def one_chromosome_rework_new(Domens, Genes, expressions, iterations) -> list:
    # Получаем домены и гены
    gene_lengths = [end - start for start, end in Genes]
    mean_gene_len = sum(gene_lengths) / len(gene_lengths) if gene_lengths else 0
    
    mu1, d1 = count_statvalues_horis(map_genes_to_domains(Domens, Genes, expressions))

    mu_list, d_list, matrix, m_v, d_v = [], [], [], [], []
    time1 = time.time()
    for _ in range(iterations):  
        shift_value = int(np.random.rand() * mean_gene_len)
        shift_value = shift_value*np.random.choice([-1, 1])
        new_arr_of_domens = rotate_genes(Domens, Genes, expressions, shift_value)
        matrix.append(new_arr_of_domens)
        
    time2 = time.time()
    for new_arr_of_domens in matrix:
        mu2, d2 = count_statvalues_horis(new_arr_of_domens)
        if mu2 is not None:
            mu_list.append(mu2)
        if d2 is not None:
            d_list.append(d2)
    mu_new = np.mean(mu_list) if mu_list else 0
    d_new = np.mean(d_list) if d_list else 0
    z = z_criterium(mu1, mu_new, d1, d_new, len(Genes)) - 1
    
    for j in range(len(matrix[0])):
        column = [row[j] for row in matrix]
        mu2, d2 = count_statvalues_horis(new_arr_of_domens)
        if mu2 is not None:
            m_v.append(mu2)
            d_v.append(d2)
        mu_vert = np.mean(m_v) if m_v else 0
        d_vert = np.mean(d_v) if d_v else 0
        
    result = [z, mu1, d1, mu_new, d_new, mu_vert, d_vert]
    return result
