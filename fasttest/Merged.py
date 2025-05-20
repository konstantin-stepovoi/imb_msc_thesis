import random
import numpy as np
import os
import multiprocessing

def process_bed_files(file1, file2):
    # Открываем файлы и читаем строки
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    # Создаем два пустых списка
    list1 = []
    list2 = []

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
                        continue  # Пропускаем строки с некорректными числами

    # Обрабатываем оба файла
    process_lines(lines1, list1)
    process_lines(lines2, list2)

    return list1, list2


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

def count_statvalues(arr_of_domens: np.ndarray) -> tuple:
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

def one_chromosome_rework(Domens, Genes, expressions) -> None:
    # Получаем домены и гены
    gene_lengths = [end - start for start, end in Genes]
    mean_gene_len = sum(gene_lengths) / len(gene_lengths) if gene_lengths else 0
    mu1, d1 = count_statvalues(map_genes_to_domains(Domens, Genes, expressions))

    mu_list, d_list = [], []
    for _ in range(5000):  
        shift_value = int(np.random.rand() * mean_gene_len)
        new_arr_of_domens = rotate_genes(Domens, Genes, expressions, shift_value)
        mu2, d2 = count_statvalues(new_arr_of_domens)
        if mu2 is not None:
            mu_list.append(mu2)
        if d2 is not None:
            d_list.append(d2)

    mu_new = np.mean(mu_list) if mu_list else 0
    d_new = np.mean(d_list) if d_list else 0
    return z_criterium(mu1, mu_new, d1, d_new, len(Genes))
    

def many_chromosome_rework(Domens, Genes, expressions):
    z_criterium_results = []
    for _ in range(1000):
        z_criterium_results.append(one_chromosome_rework(Domens, Genes, expressions))
    return z_criterium_results



def process_files(pair_index):
    # Генерация имен файлов на основе индекса пары
    gene_file = f"Genes/randomG.{pair_index:05d}.bed"
    domain_file = f"Domains/randomD.{pair_index:05d}.bed"

    # В случае ошибки при открытии файлов (например, если файл не существует)
    try:
        # Получаем Domens и Genes с помощью функции из JFL.py
        Domens, Genes = process_bed_files(domain_file, gene_file)
        expressions = [random.randint(0, 55) for _ in range(len(Genes))]

        # Обрабатываем хромосому и сохраняем результат
        results = many_chromosome_rework(Domens, Genes, expressions) 

        # Сохраняем результат в файл
        with open(f"results/{pair_index}.txt", "w") as result_file:
            for result in results:
                result_file.write(f"{result}\n")
        print(f'Processed file pair {pair_index}')

    except Exception as e:
        print(f"Error processing pair {pair_index}: {e}")

def main():
    # Папки, содержащие файлы с генами и доменами
    genes_dir = "Genes"
    domains_dir = "Domains"

    # Получаем все номера файлов (1-1000)
    files_count = len([name for name in os.listdir(genes_dir) if name.endswith('.bed')])

    # Создаем папку для результатов, если ее нет
    if not os.path.exists("results"):
        os.makedirs("results")

    # Используем multiprocessing для параллельной обработки
    # Получаем метод запуска multiprocessing
    method = multiprocessing.get_start_method()

    # Обработка всех пар файлов в параллели
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(process_files, range(1, files_count + 1))  # Обрабатываем все файлы

if __name__ == "__main__":
    multiprocessing.freeze_support()  # Поддержка для заморозки в Windows
    main()
