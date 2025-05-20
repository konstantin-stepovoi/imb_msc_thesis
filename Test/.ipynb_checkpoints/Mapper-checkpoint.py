import numpy as np

def map_genes_to_domains_3(Domens, Genes, expressions):
    arr_of_domens = [None] * len(Domens)
    gene_idx = 0

    for dom_idx, (dom_start, dom_end) in enumerate(Domens):
        domain_expressions = []
        while gene_idx < len(Genes) and Genes[gene_idx][0] < dom_start:
            gene_idx += 1
        while gene_idx < len(Genes) and Genes[gene_idx][0] <= dom_end:
            domain_expressions.append(expressions[gene_idx])
            gene_idx += 1
        
        arr_of_domens[dom_idx] = domain_expressions if domain_expressions else None

    return np.array(arr_of_domens, dtype=object)

def map_genes_to_domains_4(Domens, Genes, expressions):
    arr_of_domens = [None] * len(Domens)
    gene_idx = 0

    for dom_idx, (dom_start, dom_end) in enumerate(Domens):
        domain_expressions = []
        
        while gene_idx < len(Genes) and Genes[gene_idx][1] < dom_start:
            gene_idx += 1
        while gene_idx < len(Genes) and Genes[gene_idx][1] <= dom_end:
            domain_expressions.append(expressions[gene_idx])
            gene_idx += 1
        arr_of_domens[dom_idx] = domain_expressions if domain_expressions else None

    return np.array(arr_of_domens, dtype=object)


def map_genes_to_domains_1(Domens, Genes, expressions):
    arr_of_domens = [None] * len(Domens)
    gene_index = 0
    for domain_index, (domain_start, domain_end) in enumerate(Domens):
        expressions_for_domain = []
        while gene_index < len(Genes) and Genes[gene_index][1] < domain_start:
            gene_index += 1
        temp_index = gene_index
        while temp_index < len(Genes) and Genes[temp_index][0] < domain_end:
            gene_start, gene_end = Genes[temp_index]
            mid_point = (gene_start + gene_end) / 2
            if domain_start <= mid_point <= domain_end:
                expressions_for_domain.append(expressions[temp_index])
            temp_index += 1
        if expressions_for_domain:
            arr_of_domens[domain_index] = expressions_for_domain

    # Возвращаем результат
    return np.array(arr_of_domens, dtype=object)


def map_genes_to_domains_2(Domens, Genes, expressions):
    arr_of_domens = [None] * len(Domens)
    domain_index = 0
    gene_index = 0
    
    while domain_index < len(Domens) and gene_index < len(Genes):
        domain_start, domain_end = Domens[domain_index]
        gene_start, gene_end = Genes[gene_index]
        gene_length = gene_end - gene_start
        if gene_end < domain_start:
            gene_index += 1
        elif gene_start > domain_end:
            domain_index += 1
        else:
            overlap_start = max(gene_start, domain_start)
            overlap_end = min(gene_end, domain_end)
            overlap_length = overlap_end - overlap_start
            expression_contribution = expressions[gene_index] * (overlap_length / gene_length)

            if arr_of_domens[domain_index] is None:
                arr_of_domens[domain_index] = [expression_contribution]
            else:
                arr_of_domens[domain_index].append(expression_contribution)

            if gene_end <= domain_end:
                gene_index += 1
            else:
                domain_index += 1
    return np.array(arr_of_domens, dtype=object)
