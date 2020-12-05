import numpy as np
import pandas as pd
from itertools import combinations
import xlrd
import time


all_disease = set()
all_category = set()
disease_category = dict()
disease_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    code = each[0]
    category = each[2]
    all_disease.add(code)
    all_category.add(category)
    disease_category[code] = category
    for each1 in each[0].split(';'):
        disease_merged[each1] = each[0]

genetic_disease = set()
file = xlrd.open_workbook('../phenotype/geneAtlas_EMR.xlsx')
sheet = file.sheets()[0]
nrows = sheet.nrows
for rownum in range(1, nrows):
    row = sheet.row_values(rownum)
    icd10 = str(row[0]).split('_')[-1]
    if icd10 in disease_merged:
        genetic_disease.add(disease_merged[icd10])

df = pd.read_csv('../phenotype/comorbidity_filter.csv')
comorbidity_list = df.values.tolist()

for type in ['snp', 'gene', 'ppi', 'pathway', 'rg']:
    print(type)
    if type == 'snp':
        path1 = '../overlap/comorbidity_snp.txt'
        path2 = '../overlap/permutation/snp.csv'
    if type == 'gene':
        path1 = '../overlap/comorbidity_gene.txt'
        path2 = '../overlap/permutation/gene.csv'
    if type == 'ppi':
        path1 = '../overlap/comorbidity_ppi.txt'
        path2 = '../overlap/permutation/ppi.csv'
    if type == 'pathway':
        path1 = '../overlap/comorbidity_pathway.txt'
        path2 = '../overlap/permutation/pathway.csv'
    if type == 'rg':
        path1 = '../overlap/comorbidity_rg.txt'
        path2 = '../overlap/permutation/rg.csv'

    comorbidity_interpreted = set()
    df = pd.read_table(path1)
    list1 = df.values.tolist()
    for each in list1:
        comorbidity_interpreted.add((each[0], each[1]))

    list1 = list()
    for each in comorbidity_list:
        code1 = each[0]
        code2 = each[1]
        category1 = each[-2]
        category2 = each[-1]
        if code1 not in genetic_disease:
            continue
        if code2 not in genetic_disease:
            continue
        if (code1, code2) in comorbidity_interpreted:
            list1.append([code1, code2, category1, category2, True])
        else:
            list1.append([code1, code2, category1, category2, False])

    df = pd.DataFrame(list1, columns=['disease1', 'disease2', 'category1', 'category2', 'overlap'])

    category_order = list(all_category)
    category_order.sort()
    dict1 = dict()
    i = 0
    for each in category_order:
        dict1[each] = i
        i += 1

    category_combination = set()
    for each in category_order:
        category_combination.add((each, each))
    for each in combinations(category_order, 2):
        category_combination.add((each[0], each[1]))

    permutation_result = dict()
    for each in category_combination:
        permutation_result[each] = list()

    i = 0
    t1 = time.time()
    while i < 100001:
        print(i)
        for each in category_combination:
            category1 = each[0]
            category2 = each[1]
            df1 = df[(df['overlap']==True) & (((df['category1']==category1) & (df['category2']==category2)) | ((df['category1']==category2) & (df['category2']==category1)))]
            a = df1.shape[0]
            permutation_result[each].append(a)
        j = 0
        while j < 5:
            df['overlap'] = np.random.permutation(df['overlap'])
            j += 1
        i += 1
    print(time.time() - t1)
    list3 = list()
    list4 = list()
    for each in category_combination:
        list3.append(each[0] + '-' + each[1])
        list4.append(permutation_result[each])

    df = pd.DataFrame(list4)
    df = df.T
    df.columns = list3
    df.to_csv(path2, index=False)