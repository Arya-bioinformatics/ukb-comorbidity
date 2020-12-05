import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats import multitest

for flag in ['with', 'without']:
    print('------------\n' + flag + '\n------------')
    if flag == 'with':
        path1 = '../overlap/comorbidity_gene.txt'
        path2 = 'a.csv'
    if flag == 'without':
        path1 = '../overlap/comorbidity_gene_rmhla.txt'
        path2 = 'a.csv'

    gene_comor_count = dict()
    comorbidity_gene = dict()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            set1 = set(str2[4].split(';'))
            comorbidity_gene[(code1, code2)] = set1
            for each in set1:
                if each not in gene_comor_count:
                    gene_comor_count[each] = 0
                gene_comor_count[each] += 1
        infile.close()

    list1 = list()
    for each in gene_comor_count:
        count = gene_comor_count[each]
        list1.append([each, count])

    df = pd.DataFrame(list1, columns=['gene', 'count'])
    df = df.sort_values(by='count', ascending=False)

    top10_gene = set()
    list1 = df.values.tolist()
    for each in list1:
        gene = each[0]
        top10_gene.add(gene)
        if len(top10_gene) == 10:
            break

    set1 = set()
    for each in comorbidity_gene:
        set2 = comorbidity_gene[each]
        if len(set2 & top10_gene) != 0:
            set1.add(each)
    print(set1)

    print('------------ total comorbidity with shared genes -----------')
    print(len(comorbidity_gene))
    print('------------ total comorbidity with top 10 genes -----------')
    print(len(set1))
    print('----------------------- ratio ------------------------------')
    print(len(set1)/len(comorbidity_gene))

    bp_gene = dict()
    all_bp_gene = set()
    with open('../overlap/GO_geneset/download/c5.bp.v7.0.entrez.gmt', 'r') as infile:
        for line in infile:
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            bp = str2[0]
            set1 = set(str2[2:])
            if '' in set1:
                set1.remove('')
            if len(set1) == 0:
                continue
            if len(set1) > 200:
                continue
            bp_gene[bp] = set1
            all_bp_gene |= set1
        infile.close()

    list1 = list()
    for each in bp_gene:
        this_bp = bp_gene[each]
        a = len(top10_gene & this_bp)
        b = len(this_bp) - a
        c = len(top10_gene & all_bp_gene) - a
        d = len(all_bp_gene) - a - b - c
        [odds, p] = fisher_exact([[a, b], [c, d]], alternative='greater')
        list1.append([each, a, odds, p])

    df = pd.DataFrame(list1, columns=['bp', 'overlap gene', 'odds', 'p'])
    df = df[df['odds']>1]
    q_list = multitest.fdrcorrection(df['p'])[1]
    df['q'] = q_list
    df1 = df.loc[df['q'] < 0.05]
    print(df1.shape[0])
    df1.to_csv(path2, index=False)