import random
import numpy as np

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
def map_genes_to_domains(Domens, Genes, expressions):
    # Читаем режим обработки из файла
    try:
        with open("map_param.txt", "r") as f:
            mode = int(f.read().strip())
    except (FileNotFoundError, ValueError):
        raise ValueError("Ошибка: файл map_param.txt не найден или содержит некорректные данные")

    if mode not in {1, 2, 3, 4}:
        raise ValueError("Некорректный режим в map_param.txt. Ожидалось 1, 2, 3 или 4.")
    
    # Создаем пустой массив для доменов
    arr_of_domens = np.empty(len(Domens), dtype=object)
    
    # Проходим по каждому домену
    for domain_index, (dom_start, dom_end) in enumerate(Domens):
        this_domen = []

        for gene_index, (gene_start, gene_end) in enumerate(Genes):
            gene_length = gene_end - gene_start
            expression = expressions[gene_index]

            # Проверяем пересечение гена с доменом
            if gene_end <= dom_start or gene_start >= dom_end:
                continue  # Ген вне домена

            if mode == 1:
                # 1. Отнести ген целиком к тому домену, в котором находится большая часть гена
                overlap_left = max(dom_start, gene_start)
                overlap_right = min(dom_end, gene_end)
                overlap_length = overlap_right - overlap_left
                if overlap_length >= gene_length / 2:
                    this_domen.append(expression)
            
            elif mode == 2:
                # 2. Отнести экспрессию гена пропорционально длине, попадающей в домен
                overlap_left = max(dom_start, gene_start)
                overlap_right = min(dom_end, gene_end)
                overlap_length = overlap_right - overlap_left
                this_domen.append(expression * (overlap_length / gene_length))
            
            elif mode == 3:
                # 3. Отнести экспрессию гена целиком к точке старта
                if gene_start >= dom_start and gene_start < dom_end:
                    this_domen.append(expression)
            
            elif mode == 4:
                # 4. Отнести экспрессию гена целиком к точке конца транскрипции
                if gene_end > dom_start and gene_end <= dom_end:
                    this_domen.append(expression)
        
        arr_of_domens[domain_index] = np.array(this_domen)

    return arr_of_domens


def rotate_genes(Domens, Genes, expressions, shift_value=None):
    # Определяем случайный сдвиг и направление, если shift_value не задан
    if shift_value is None:
        max_shift = max(gene[1] for gene in Genes)
        shift_value = np.random.randint(1, max_shift) * np.random.choice([-1, 1])

    # Длина хромосомы
    chrom_length = max(gene[1] for gene in Genes)

    # Сдвигаем координаты всех генов
    new_genes = []
    for start, end in Genes:
        new_start = (start + shift_value) % chrom_length
        new_end = (end + shift_value) % chrom_length
        
        # Если ген "разрывается", корректно переносим его на новую позицию
        if new_end < new_start:
            gene_length = end - start  # Длина гена остаётся неизменной
            new_start = 0
            new_end = gene_length

        new_genes.append((new_start, new_end))

    # Сортируем гены после сдвига
    new_genes.sort()
    # Циклический сдвиг массива expressions
    shift_mod = shift_value % len(expressions)
    new_expressions = np.roll(expressions, shift_mod)
    # Вызываем функцию маппинга с новыми сдвинутыми данными
    arr_of_domens = map_genes_to_domains(Domens, new_genes, new_expressions)

    return arr_of_domens

def count_statvalues_horis(arr_of_domens: np.ndarray) -> tuple:
    means = []
    variances = []

    for domen in arr_of_domens:
        if domen is not None and len(domen) > 0:
            means.append(np.mean(domen))
            variances.append(np.var(domen))
    
    overall_mean = np.mean(means) if means else None
    overall_dispersion = np.mean(variances) if variances else None

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

    return arr_of_domens


def one_chromosome_rework_new(Domens, Genes, expressions, iterations) -> None:
    # Получаем домены и гены
    gene_lengths = [end - start for start, end in Genes]
    mean_gene_len = sum(gene_lengths) / len(gene_lengths) if gene_lengths else 0
    mu1, d1 = count_statvalues_horis(map_genes_to_domains(Domens, Genes, expressions))

    mu_list, d_list, matrix, m_v, d_v = [], [], [], [], []
    for _ in range(iterations):  
        shift_value = int(np.random.rand() * mean_gene_len)
        new_arr_of_domens = rotate_genes(Domens, Genes, expressions, shift_value)
        matrix.append(new_arr_of_domens)
        
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

    
iterations = [500, 1000, 5000, 10000]
def many_chromosome_rework(Domens, Genes, expressions):
    z_criterium_results = []
    for i in range(len(iterations)):  # Выполняем NUM итераций, как указано
        z_criterium_results.append(one_chromosome_rework_new(Domens, Genes, expressions, iterations[i]))
    return z_criterium_results