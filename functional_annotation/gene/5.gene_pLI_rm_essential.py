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


lethality_phenotype = set()
with open('../overlap/MGI_download/VOC_MammalianPhenotype.rpt.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        if 'lethal' in line:
            if ('embryonic' in line) | ('prenatal' in line) | ('postnatal' in line):
                lethality_phenotype.add(str2[0])
    infile.close()

gene_lethality = set()
with open('../overlap/MGI_download/MGI_PhenoGenoMP.rpt.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        if not str2[3].startswith('MP'):
            continue
        phenotype = str2[3]
        if phenotype not in lethality_phenotype:
            continue
        gene = str2[5]
        gene_lethality.add(gene)
    infile.close()

gene_valid_knockout = set()
with open('../overlap/MGI_download/MGI_PhenotypicAllele.rpt.txt', 'r') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        gene = str2[6]
        if gene not in gene_lethality:
            continue
        type = str2[4]
        if 'knockout' not in type:
            continue
        phenotype = str2[10]
        if not phenotype.startswith('MP'):
            continue
        set1 = set(phenotype.split(','))
        # 'MP:0010768' - mortality/aging	the observable characteristics related to the ability of a mammalian organism to live and age that are manifested throughout development and life span
        # this is due to hierarchical structure
        if 'MP:0010768' not in set1:
            continue
        gene_valid_knockout.add(gene)
    infile.close()

print(len(gene_valid_knockout))

essential_gene = set()
with open('../overlap/MGI_download/HMD_HumanPhenotype.rpt.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        entrezid = str2[1].strip(' ')
        mgi_id = str2[5].strip(' ')
        if mgi_id not in gene_valid_knockout:
            continue
        phenotype = str2[6].strip(' ')
        set1 = set(phenotype.split(' '))
        if 'MP:0010768' not in set1:
            continue
        essential_gene.add(entrezid)
    infile.close()

print('total essential gene:' + str(len(essential_gene)))



gene_pLI = dict()
with open('../overlap/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data_geneid_format.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        gene = str2[1]
        if gene in essential_gene:
            continue
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