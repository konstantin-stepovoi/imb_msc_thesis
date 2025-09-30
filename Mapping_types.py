CALIBRATION_WEIGHTS = [] 
GENE_COUNT = 0
        
def map_genes_to_domains_start(Domens, Genes, expressions): 
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    for dom_start, dom_end in Domens:
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            if dom_start <= gene_start < dom_end:
                this_domen.append(expressions[i])

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.ones(len(this_domen), dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens


def map_genes_to_domains_middle(Domens, Genes, expressions): 
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    for dom_start, dom_end in Domens:
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            gene_middle = (gene_start + gene_end) // 2
            if dom_start <= gene_middle < dom_end:
                this_domen.append(expressions[i])

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.ones(len(this_domen), dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens

def map_genes_to_domains_stop(Domens, Genes, expressions): 
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    for dom_start, dom_end in Domens:
        this_domen = []
        for i, (gene_start, gene_end) in enumerate(Genes):
            if dom_start < gene_end <= dom_end:
                this_domen.append(expressions[i])

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.ones(len(this_domen), dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens



def affine_mod(sig1, sig2):
    r = (sig1 - sig2) / sig1
    lower, upper = 0.076, 0.105
    scale = (upper - lower) / 8.5 
    center = (upper + lower) / 2
    z = (r - center) / scale + 0.3
    return z
    
def map_genes_to_domains_propor(Domens, Genes, expressions):
    global CALIBRATION_WEIGHTS, GENE_COUNT, DOMAIN_LENGTHS

    GENE_COUNT = len(Genes)
    valid_domens = []
    valid_weights = []

    gene_index = 0
    for dom_start, dom_end in Domens:
        this_domen = []
        this_weights = []

        while gene_index < len(Genes) and Genes[gene_index][1] < dom_start:
            gene_index += 1
        temp_index = gene_index

        while temp_index < len(Genes) and Genes[temp_index][0] < dom_end:
            gene_start, gene_end = Genes[temp_index]
            gene_length = gene_end - gene_start
            expression = expressions[temp_index]

            weight = 0.0
            if dom_start <= gene_start and gene_end <= dom_end:
                weight = 1.0
            elif gene_start <= dom_start and dom_end <= gene_end:
                weight = (dom_end - dom_start) / gene_length
            elif gene_start < dom_start < gene_end <= dom_end:
                weight = (gene_end - dom_start) / gene_length
            elif dom_start <= gene_start < dom_end < gene_end:
                weight = (dom_end - gene_start) / gene_length

            if weight > 0:
                this_domen.append(expression * weight)
                this_weights.append(weight)

            temp_index += 1

        if len(this_domen) > 1:
            valid_domens.append(np.array(this_domen, dtype=float))
            valid_weights.append(np.array(this_weights, dtype=float))

    arr_of_domens = np.array(valid_domens, dtype=object)
    CALIBRATION_WEIGHTS = np.array(valid_weights, dtype=object)
    DOMAIN_LENGTHS = np.array([len(d) for d in arr_of_domens], dtype=int)

    return arr_of_domens