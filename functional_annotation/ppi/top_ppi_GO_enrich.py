import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats import multitest
from collections import Counter


for flag in ['with', 'without']:
    print('------------\n' + flag + '\n------------')
    if flag == 'with':
        path1 = '../overlap/comorbidity_ppi.txt'
        path2 = 'a.csv'
    if flag == 'without':
        path1 = '../overlap/comorbidity_ppi_rmhla.txt'
        path2 = 'a1.csv'

    ppi_comor_count = dict()
    ppi_comorbidity = dict()
    all_ppi_interpreted = 0
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            set1 = set(str2[4].split(';'))
            all_ppi_interpreted += 1
            for each in set1:
                if each not in ppi_comor_count:
                    ppi_comor_count[each] = 0
                    ppi_comorbidity[each] = set()
                ppi_comor_count[each] += 1
                ppi_comorbidity[each].add((code1, code2))
        infile.close()

    df = pd.DataFrame.from_dict(ppi_comor_count, orient='index', columns=['count'])
    df['ppi'] = df.index
    df = df.sort_values(by='count', ascending=False)
    list1 = df['ppi'].values.tolist()
    i = 1
    top_ppi_gene = set()
    top_ppi = set()
    for each in list1:
        top_ppi.add(each)
        [gene1, gene2] = each.split('~')
        top_ppi_gene.add(gene1)
        top_ppi_gene.add(gene2)
        if i == 10:
            break
        i += 1

    bp_gene = dict()
    all_bp_gene = set()
    i = 0
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
                i += 1
                continue
            bp_gene[bp] = set1
            all_bp_gene |= set1
        infile.close()

    list1 = list()
    for each in bp_gene:
        this_bp = bp_gene[each]
        a = len(top_ppi_gene & this_bp)
        b = len(this_bp) - a
        c = len(top_ppi_gene & all_bp_gene) - a
        d = len(all_bp_gene) - a - b - c
        [odds, p] = fisher_exact([[a, b], [c, d]], alternative='greater')
        list1.append([each, a, odds, p])

    df = pd.DataFrame(list1, columns=['bp', 'overlap gene', 'odds', 'p'])
    df = df[df['odds'] > 1]
    q_list = multitest.fdrcorrection(df['p'])[1]
    df['q'] = q_list
    df1 = df.loc[df['q'] < 0.05]
    print(df1.shape[0])
    df1.to_csv(path2, index=False)

    set2 = set()
    for each in top_ppi:
        set1 = ppi_comorbidity[each]
        set2 |= set1
    list1 = list()
    for each in set2:
        list1.append(each[0])
        list1.append(each[1])
    x = Counter(list1)
    print(set2)
    list1 = list()
    for each in set2:
        list1.append(each[0])
        list1.append(each[1])
    result = Counter(list1)
    print('all_ppi_interpreted comorbiidty: ' + str(all_ppi_interpreted))
    print('top ppi involved comorbidity: ' + str(len(set2)))
    print('ratio: ' + str(len(set2)/all_ppi_interpreted))
