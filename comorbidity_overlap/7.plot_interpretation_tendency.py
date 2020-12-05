import pandas as pd
import numpy as np
from statsmodels.stats import multitest
import xlrd
import os


disease_class = dict()
all_category = set()
disease_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    code = each[0]
    category = each[2]
    disease_class[code] = category
    all_category.add(category)
    for each1 in code.split(';'):
        disease_merged[each1] = each[0]

category_order = list(all_category)
category_order.sort()
category_index = dict()
i = 0
for each in category_order:
    category_index[each] = i
    i += 1

genetic_disease = set()
file = xlrd.open_workbook('../phenotype/geneAtlas_EMR.xlsx')
sheet = file.sheets()[0]
nrows = sheet.nrows
for rownum in range(1, nrows):
    row = sheet.row_values(rownum)
    icd10 = str(row[0]).split('_')[-1]
    if icd10 in disease_merged:
        icd10 = disease_merged[icd10]
        genetic_disease.add(icd10)
print(len(genetic_disease))

genetic_comorbidity = list()
df = pd.read_csv('../phenotype/comorbidity_filter.csv')
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    if (code1 in genetic_disease) & (code2 in genetic_disease):
        c1 = disease_class[code1]
        c2 = disease_class[code2]
        genetic_comorbidity.append([code1, code2, c1, c2])

genetic_df = pd.DataFrame(genetic_comorbidity, columns=['code1', 'code2', 'category1', 'category2'])

significant_set = set()
category_category_sig_combination = dict()
for each in ['snp', 'gene', 'ppi', 'pathway', 'rg']:
    print(each)
    if each == 'snp':
        path1 = '../overlap/permutation/snp.csv'
        path2 = 'snp.pdf'
        # continue
    elif each == 'gene':
        path1 = '../overlap/permutation/gene.csv'
        path2 = 'gene.pdf'
        # continue
    elif each == 'ppi':
        path1 = '../overlap/permutation/ppi.csv'
        path2 = 'ppi.pdf'
        # continue
    elif each == 'pathway':
        path1 = '../overlap/permutation/pathway.csv'
        path2 = 'pathway.pdf'
        # continue
    elif each == 'rg':
        path1 = '../overlap/permutation/rg.csv'
        path2 = 'rg.pdf'
        # continue

    df = pd.read_csv(path1)
    list1 = df.values.tolist()
    category_combination = df.columns

    list1 = list()
    for each1 in category_combination:
        [category1, category2] = each1.split('-')
        a = df[each1][0]
        cond = (((genetic_df['category1'] == category1) & (genetic_df['category2'] == category2)) | (
                (genetic_df['category2'] == category1) & (genetic_df['category1'] == category2)))
        df1 = genetic_df[cond]
        b = df1.shape[0]
        df1 = df[df[each1] >= a]
        c = df1.shape[0]
        if b == 0:
            list1.append([each1, a, 0, float(c) / df.shape[0]])
        else:
            list1.append([each1, a, float(a) / b, float(c) / df.shape[0]])

    df = pd.DataFrame(list1, columns=['category_category', 'interperted count', 'interpreted ratio', 'p'])
    q_list = multitest.fdrcorrection(df['p'])[1]
    df['q'] = q_list
    significant_set |= set(df.loc[df['q']<0.05, 'category_category'])
    for each1 in set(df.loc[df['q']<0.05, 'category_category']):
        if each1 not in category_category_sig_combination:
            category_category_sig_combination[each1] = list()
        if each == 'snp':
            category_category_sig_combination[each1].append('1')
        elif each == 'gene':
            category_category_sig_combination[each1].append('2')
        elif each == 'ppi':
            category_category_sig_combination[each1].append('3')
        elif each == 'pathway':
            category_category_sig_combination[each1].append('4')
        elif each == 'rg':
            category_category_sig_combination[each1].append('5')


    sig_dict = dict()
    ratio_data = np.zeros((len(category_order), len(category_order)))
    q_data = np.zeros((len(category_order), len(category_order)))

    list1 = df.values.tolist()
    for each1 in list1:
        [category1, category2] = each1[0].split('-')
        count = each1[1]
        ratio = each1[2]
        q = each1[-1]
        i = category_index[category1]
        j = category_index[category2]
        ratio_data[i, j] = ratio
        ratio_data[j, i] = ratio
        q_data[i, j] = q
        q_data[j, i] = q

    df = pd.DataFrame(ratio_data)
    df.columns = category_order
    df.index = category_order
    df.to_csv('ratio.csv')

    df = pd.DataFrame(q_data)
    df.columns = category_order
    df.index = category_order
    df.to_csv('q.csv')

    os.system('Rscript 7.plot_corr.R' + ' ' + path2)

m = 0
n = 0
for each in significant_set:
    c1, c2 = each.split('-')
    if c1 == c2:
        m += 1
    else:
        n += 1
print('significant interpreted intra-category: ' + str(m))
print('significant interpreted inter-category: ' + str(n))


comb_data = np.zeros((len(category_order), len(category_order)))
for each in category_category_sig_combination:
    [category1, category2] = each.split('-')
    comb = int(''.join(category_category_sig_combination[each]))
    i = category_index[category1]
    j = category_index[category2]
    comb_data[i, j] = comb
    comb_data[j, i] = comb
df = pd.DataFrame(comb_data)
df.columns = category_order
df.index = category_order
df.to_csv('c.csv')