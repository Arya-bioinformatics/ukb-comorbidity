import pandas as pd


all_gene = set()
with open('../genome/grch37_download/genome_assemblies/gene_assembly_grch37.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\n')
        str2 = str1.split('\t')
        if str2[5] == 'pseudo':
            continue
        all_gene.add(str2[0])
    infile.close()

ppi_dict = dict()
with open('../genome/PPI_biogrid/download/BIOGRID-ALL-3.5.176.tab2.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        gene1 = str2[1]
        gene2 = str2[2]
        orgnism_A = str2[15]
        orgnism_B = str2[16]
        experimental_system = str2[11]
        if (orgnism_A != '9606') | (orgnism_B != '9606') | (gene1 not in all_gene) | (gene2 not in all_gene):
            continue
        if gene1 not in ppi_dict:
            ppi_dict[gene1] = set()
        ppi_dict[gene1].add(gene2)
        if gene2 not in ppi_dict:
            ppi_dict[gene2] = set()
        ppi_dict[gene2].add(gene1)
    infile.close()

for flag in ['with', 'without']:
    print(' ----------------- ' + flag + '----------------- ')
    if flag == 'with':
        path1 = '../genome/disease_gene.txt'
        path2 = '../genome/disease_ppi.txt'
    elif flag == 'without':
        path1 = '../genome/disease_gene_rmhla.txt'
        path2 = '../genome/disease_ppi_rmhla.txt'

    disease_gene = dict()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_gene[str2[0]] = set(str2[1].split(';'))
        infile.close()

    disease_ppi = dict()
    for each in disease_gene:
        gene_set = disease_gene[each]
        ppi_set = set()
        for each1 in gene_set:
            if each1 in ppi_dict:
                set1 = ppi_dict[each1]
                for each2 in set1:
                    if each1 < each2:
                        ppi_set.add(each1 + '~' + each2)
                    else:
                        ppi_set.add(each2 + '~' + each1)
        if len(ppi_set) == 0:
            continue
        disease_ppi[each] = ppi_set

    print('total disease with ppi: ' + str(len(disease_ppi)))

    df = pd.read_csv('../phenotype/comorbidity_filter.csv')
    comorbidity_ppi = dict()
    list1 = df.values.tolist()
    for each in list1:
        code1, code2 = each[:2]
        if (code1 in disease_ppi) & (code2 in disease_ppi):
            set1 = disease_gene[code1]
            set2 = disease_gene[code2]
            set3 = disease_ppi[code1] & disease_ppi[code2]
            set4 = set()
            for each1 in set3:
                gene1, gene2 = each1.split('~')
                if (len(set([gene1, gene2]) & set1) == 2) | (len(set([gene1, gene2]) & set2) == 2):
                    continue
                if ((gene1 in set1) & (gene2 in set2)) | ((gene2 in set1) & (gene1 in set2)):
                    set4.add(each1)
            if len(set4) != 0:
                comorbidity_ppi[(code1, code2)] = set4

    print('total comorbidity with shared ppi: ' + str(len(comorbidity_ppi)))