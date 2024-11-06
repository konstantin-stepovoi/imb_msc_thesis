from dataclasses import dataclass

@dataclass
class Config:
    '''этот класс содержит только параметры моделей, которые мы используем
    и параметры вычислительного ядра
    его можно менять под задачи, а старые значения не потеряются потому что они 
    бэкаются в outputах
    '''
    #глобольные параметры программы:
    rotations_per_sample: int = 10000
    memorytype: str = 'file' #есть ещё 'buffer'
    domentype: str = 'uni' #есть ещё 'lognorm', 'uni'
    genetype = 'exp' #есть ещё 'exp', 'uni'
    chromosome_length = 800000000
    
    
    # Параметры модели распределения длин генов
    y0: float = 2.18
    x0: float = -356.69
    A1: float = 15019.03
    t1: float = 93.14
    A2: float = 16085.19
    t2: float = 93.14
    A3: float = 67.25
    t3: float = 2165.97
    Min_gene_len: int = 50
    Max_gene_len: int = 66000
    mean_gene_length: int = 34000

    # Параметры модели генов
    number_of_genes: int = 5000
    minimal_expression: int = 0
    maximal_expression: int = 50
    
    #Параметры модели доменов
    domen_length: int = 110000 
    min_dom_length: int = 10000
    max_dom_length: int = 1000000
    
