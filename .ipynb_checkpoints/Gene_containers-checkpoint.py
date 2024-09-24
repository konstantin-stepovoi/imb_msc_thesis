import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional
'''
этот файл просто определяет что такое ген и что такое хромосома.
сюда будем добавлять новые датаклассы потом.
этот файл используется Modeller'ом, чтобы генерировать новые данные, и Parser'ом,
чтобы приводить реальные данные к удобной форме.
'''
@dataclass
class Gene:
    length: np.int32
    expression: np.uint8
    start_point: np.int32
    end_point: np.int32
    name: Optional[str] = field(default=None)

    def __post_init__(self):
        if self.length != self.end_point - self.start_point:
            raise ValueError("Length does not match the difference between start and end points.")
class Chromosome(np.ndarray):
    def __new__(cls, genes: List[Gene] = None):
        if genes is None:
            genes = []
        obj = np.asarray(genes, dtype=object).view(cls)
        return obj
    
    def add_gene(self, gene: Gene) -> 'Chromosome':
        new_size = len(self) + 1
        new_array = np.empty(new_size, dtype=object)
        new_array[:-1] = self
        new_array[-1] = gene
        return new_array.view(Chromosome)
    
    def index(self, search_value: int, by: str = 'start_point') -> int:
        for i, gene in enumerate(self):
            if hasattr(gene, by) and getattr(gene, by) == search_value:
                return i
        raise ValueError(f"{search_value} not found in Chromosome by {by}")

    def rotate(self, shift_value: int) -> 'Chromosome':
        chromosome_length = self[-1].end_point
        for gene in self:
            gene.start_point += shift_value
            if gene.start_point > chromosome_length:
                gene.start_point -= chromosome_length
            gene.end_point += shift_value
            if gene.end_point > chromosome_length:
                gene.end_point -= chromosome_length

        indices = np.argsort([gene.start_point for gene in self])
        self[:] = self[indices]
        return self