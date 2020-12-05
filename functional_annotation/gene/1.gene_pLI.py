import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import seaborn as sns



all_gene= set()
with open('../genome/grch37_download/genome_assemblies/gene_assembly_grch37.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\n')
        str2 = str1.split('\t')
        gene_id = str2[0]
        if type == 'pseudo':
            continue
        all_gene.add(gene_id)
    infile.close()

gene_pLI = dict()
with open('../overlap/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data_geneid_format.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        gene = str2[1]
        pLI = float(str2[19])
        gene_pLI[gene] = pLI
    infile.close()

for flag in ['with', 'without']:
    print('\n' + flag + '\n')
    if flag == 'with':
        path1 = '../genome/disease_gene.txt'
        path2 = '../overlap/comorbidity_gene.txt'
        path3 = 'a.pdf'
        # continue
    if flag == 'without':
        path1 = '../genome/disease_gene_rmhla.txt'
        path2 = '../overlap/comorbidity_gene_rmhla.txt'
        path3 = 'a1.pdf'
        # continue

    disease_gene = set()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_gene = disease_gene | set(str2[1].split(';'))
        infile.close()

    comorbidity_gene = set()
    with open(path2, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            comorbidity_gene = comorbidity_gene | set(str2[4].split(';'))
        infile.close()

    nonassociated_gene_pLI = list()
    only_disease_gene_pLI = list()
    comorbidity_gene_pLI = list()
    for each in gene_pLI:
        pLI = gene_pLI[each]
        if each in comorbidity_gene:
            comorbidity_gene_pLI.append(pLI)
        elif each in disease_gene:
            only_disease_gene_pLI.append(pLI)
        elif each in all_gene:
            nonassociated_gene_pLI.append(pLI)

    [u1, p1] = ttest_ind(comorbidity_gene_pLI, nonassociated_gene_pLI)
    [u2, p2] = ttest_ind(comorbidity_gene_pLI, only_disease_gene_pLI)
    [u3, p3] = ttest_ind(only_disease_gene_pLI, nonassociated_gene_pLI)

    print([u1, p1])
    print([u2, p2])
    print([u3, p3])

    list1 = list()
    for each in nonassociated_gene_pLI:
        list1.append(['Non-associated', each])
    for each in only_disease_gene_pLI:
        list1.append(['Only disease', each])
    for each in comorbidity_gene_pLI:
        list1.append(['Comorbidity', each])

    df = pd.DataFrame(list1, columns=['group', 'pLI'])
    b = sns.violinplot(x="group", y="pLI", data=df)
    b.set_xlabel("", fontsize=1)
    b.set_ylabel("pLI", fontsize=16)
    b.tick_params(labelsize=16)
    b.set_xticklabels(['Non-disease genes', 'Other disease-genes', 'Comorbidity-genes'], rotation=45)
    plt.savefig(path3, bbox_inches='tight')
    plt.close()