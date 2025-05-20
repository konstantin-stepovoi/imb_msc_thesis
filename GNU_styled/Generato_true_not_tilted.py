#!/usr/bin/python3
import pandas as pd
import os
import random
import numpy as np

def make_new_files_for_run(genes_path, domains_path, file_number, percentage):
    unique_genes = set()
    genes = []
    percentage = percentage
    with open('genes_reference.txt', 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue  
            _, start, end, expression = parts
            start = int(start)
            end = int(end)
            expression = float(expression)
            center = (start + end) // 2
    
            if center not in unique_genes:
                genes.append((start, end, expression, center))
                unique_genes.add(center)
                
    genes.sort(key=lambda x: x[3])
    df = pd.read_excel('Domains_reference.xlsx', header=None)  # если нет заголовков
    df.columns = ['chrom', 'start', 'end']
    chr3_domains = df[df['chrom'] == 'chr3']
    domain_coords = list(zip(chr3_domains['start'], chr3_domains['end']))
    
    results = []
    
    for domain_start, domain_end in domain_coords:
        genes_in_domain = [g for g in genes if domain_start <= g[3] <= domain_end]
        count = len(genes_in_domain)
        if count > 0:
            avg_expression = sum(g[2] for g in genes_in_domain) / count
        else:
            avg_expression = 0.0
    
        results.append((domain_start, domain_end, count, avg_expression))
    
    domains_path = domains_path
    genes_path = genes_path
    file_number = file_number
    
    os.makedirs(domains_path, exist_ok=True)
    os.makedirs(genes_path, exist_ok=True)
    
    chr_name = "chr3"
    bed_lines = [f"{chr_name}\t{start}\t{end}\n" for (start, end, _, _) in results]
    
    for i in range(0, file_number):
        filename = f"{domains_path}/randomD.{i:05d}.bed"
        with open(filename, "w") as f:
            f.writelines(bed_lines)
    
    for file_idx in range(0, file_number):
        filename = f"{genes_path}/randomG.{file_idx:05d}.bed"
        with open(filename, "w") as f:
            for dom_start, dom_end, gene_count, avg_expr in results:
                gene_count = int(gene_count)
                if gene_count == 0:
                    continue
    
                domain_length = dom_end - dom_start
                spacing = 10  # отступ между генами и от краёв
    
                usable_length = domain_length - spacing * (gene_count + 1)
                if usable_length <= gene_count:
                    continue
    
                gene_lengths = [usable_length // gene_count] * gene_count
                gene_lengths[-1] += usable_length % gene_count  # добиваем остаток
    
                current_start = dom_start + spacing
                for gene_len in gene_lengths:
                    gene_end = current_start + gene_len
                    # expression = avg_expr * random.uniform(1-percentage, 1+percentage)
                    expression = np.random.randint(0, 50)
                    f.write(f"{chr_name}\t{int(current_start)}\t{int(gene_end)}\t{expression:.3f}\t+\n")
                    current_start = gene_end + spacing