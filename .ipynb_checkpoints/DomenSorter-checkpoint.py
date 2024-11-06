import numpy as np
from Gene_containers import Gene, Chromosome
from parameters import Config, glob_config
config = Config()
'''
этот файл принимает хромосому и отдаёт массив доменов,
каждый из элементов которого заполнен согласно определяемыми 
нами правилами и файлу параметров
'''

def generate_lognormal_value(median: float, lower_bound: float, upper_bound: float) -> float:
    mu = np.log(median)
    z_lower = (np.log(lower_bound) - mu)
    z_upper = (np.log(upper_bound) - mu)
    sigma = (z_upper - z_lower) / (2 * 3.91993)
    return np.random.lognormal(mean=mu, sigma=sigma)



def create_empty_object_array(size: int) -> np.ndarray:
    return np.empty(size, dtype=object)

def add_array(array_of_arrays: np.ndarray, index: int, new_array: np.ndarray):
    array_of_arrays[index] = new_array

def Sort_to_domens(chromosome: Chromosome, domlen: int) -> np.ndarray:
    number_of_domens = (chromosome[-1].end_point // domlen) + 1
    arr_of_domens = create_empty_object_array(number_of_domens)
    
    gap_size = config.gapsizes[config.cur_gs_index] # Получаем текущий размер гэпа

    domlens_r = []
    domlens_l = []
    for n in range(number_of_domens):
        if config.domentype == 'uni':
            domlens_r.append((domlen + gap_size)*n + domlen)
            domlens_l.append((domlen * n) + gap_size)
        elif config.domentype == 'lognorm':
            generated_length = generate_lognormal_value(median=config.domen_length,
                                                         lower_bound=config.min_dom_length, upper_bound=config.max_dom_length)
            start_point = (domlens_r[-1] if domlens_r else 0) + gap_size  # Начало домена с учетом предыдущей правой границы и гэпа
            end_point = start_point + generated_length  # Конец домена
            domlens_l.append(start_point)
            domlens_r.append(end_point)
        else:
            raise ValueError("Invalid type: choose 'uni' or 'lognorm'")
    
    for n in range(number_of_domens):
        this_domen = []
        this_d_l = domlens_l[n]
        this_d_r = domlens_r[n]
    
        for gene in chromosome:
            if this_d_l <= gene.start_point < this_d_r and this_d_l < gene.end_point <= this_d_r:
                this_domen.append(gene.expression)
            elif gene.start_point < this_d_l < gene.end_point:
                overlap_length = gene.end_point - this_d_l
                this_domen.append(gene.expression * overlap_length / gene.length)
            elif this_d_l <= gene.start_point < this_d_r and gene.end_point > this_d_r:
                overlap_length = this_d_r - gene.start_point
                this_domen.append(gene.expression * overlap_length / gene.length)
    
        if this_domen:
            add_array(arr_of_domens, n, np.array(this_domen))
    
    return arr_of_domens

def Group_Chromosomes_by_lines(file_name: str):
    chromosomes = {}
    
    with open(file_name, 'r') as file:
        lines = file.readlines()
        
        for line in lines:
            # Разбиваем строку на составляющие
            data = line.split()
            chrom = data[0]  # Метка хромосомы (chr0, chr1 и т.д.)
            
            # Если хромосомы ещё нет в словаре, создаём её
            if chrom not in chromosomes:
                chromosomes[chrom] = []
                
            # Добавляем строку в соответствующую хромосому
            chromosomes[chrom].append(line.strip())
    return chromosomes

def Domen_Sort_f(file_name: str, chromnum: int) -> np.ndarray:
    chromosomes = Group_Chromosomes_by_lines(file_name)
    current_chromosome = chromosomes.get(f'chr{chromnum}', [])

    if not current_chromosome:
        raise ValueError(f"No data found for chromosome chr{chromnum}")

    # Создаём объект Chromosome на основе данных
    genes = []
    for line in current_chromosome:
        parts = line.split()
        if len(parts) != 6:
            raise ValueError(f"Invalid line format: {line}")
        
        # Извлекаем данные
        start_point = np.int32(parts[1])
        end_point = np.int32(parts[2])
        name = parts[3]
        expression = np.uint8(parts[4])  # Используем 4-й элемент как expression
        length = end_point - start_point
        
        # Создаем экземпляр Gene
        gene = Gene(length=length, expression=expression, start_point=start_point, end_point=end_point, name=name)
        genes.append(gene)
    
    # Создаем Chromosome из списка genes
    chromosome = Chromosome(genes)
    domenlength = config.domen_length
    arr_of_domens = Sort_to_domens(chromosome, domenlength)
    return arr_of_domens
    
    
    
