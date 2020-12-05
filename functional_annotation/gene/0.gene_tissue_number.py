import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu


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


df = pd.read_table('../GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_tissue_mean_geneid_format.txt')
df['geneid'] = df['geneid'].astype(str)
all_tissue_gene = set(df['geneid']) & all_gene
tissue_expressed_genes = dict()
all_tissue = list(df.columns)[1:]
for each in all_tissue:
    df1 = df.loc[df[each] > 1, 'geneid']
    set1 = set(df1)
    if len(set1) == 0:
        continue
    tissue_expressed_genes[each] = set1 & all_gene


all_gene_tissue_number = dict()
for each in all_gene:
    i = 0
    for each1 in tissue_expressed_genes:
        set1 = tissue_expressed_genes[each1]
        if each in set1:
            i += 1
    if i == 0:
        continue
    all_gene_tissue_number[each] = i


for flag in ['with', 'without']:
    print('\n' + flag + '\n')
    if flag == 'with':
        path1 = '../genome/disease_gene.txt'
        path2 = '../overlap/comorbidity_gene.txt'
        path3 = 'a.csv'
        # continue
    if flag == 'without':
        path1 = '../genome/disease_gene_rmhla.txt'
        path2 = '../overlap/comorbidity_gene_rmhla.txt'
        path3 = 'a1.csv'
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

    nonassociated_gene_tissues = list()
    only_disease_gene_tissues = list()
    comorbidity_gene_tissues = list()
    for each in all_gene_tissue_number:
        num = all_gene_tissue_number[each]
        if each in comorbidity_gene:
            comorbidity_gene_tissues.append(num)
        elif each in disease_gene:
            only_disease_gene_tissues.append(num)
        else:
            nonassociated_gene_tissues.append(num)

    [u1, p1] = mannwhitneyu(comorbidity_gene_tissues, nonassociated_gene_tissues)
    [u2, p2] = mannwhitneyu(comorbidity_gene_tissues, only_disease_gene_tissues)
    [u3, p3] = mannwhitneyu(only_disease_gene_tissues, nonassociated_gene_tissues)

    print([u1, p1])
    print([u2, p2])
    print([u3, p3])

    s1 = pd.Series(nonassociated_gene_tissues)
    s2 = pd.Series(only_disease_gene_tissues)
    s3 = pd.Series(comorbidity_gene_tissues)
    df = pd.DataFrame()
    df['Non-associated'] = s1
    df['Only disease'] = s2
    df['Comorbidity'] = s3

    r_list = [0, 1, 5, 42, 53]
    with open(path3, 'w+') as outfile:
        outfile.write('Tissue number,Non-associated,Only disease,Comorbidity\n')
        for i in range(0, len(r_list) - 1):
            df1 = df.loc[(df['Non-associated'] > r_list[i]) & (df['Non-associated'] <= r_list[i + 1])]
            a = df1.shape[0] / float(len(nonassociated_gene_tissues))
            df1 = df.loc[(df['Only disease'] > r_list[i]) & (df['Only disease'] <= r_list[i + 1])]
            b = df1.shape[0] / float(len(only_disease_gene_tissues))
            df1 = df.loc[(df['Comorbidity'] > r_list[i]) & (df['Comorbidity'] <= r_list[i + 1])]
            c = df1.shape[0] / float(len(comorbidity_gene_tissues))
            outfile.write(
                str(r_list[i] + 1) + '-' + str(r_list[i + 1]) + ',' + str(a) + ',' + str(b) + ',' + str(c) + '\n')
    outfile.close()