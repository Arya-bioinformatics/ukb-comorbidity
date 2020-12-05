import pandas as pd
from itertools import combinations
from statsmodels.stats import multitest
from itertools import product
import numpy as np
import rpy2.robjects as robjects
from matplotlib import pyplot as plt


disease_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    for each1 in each[0].split(';'):
        disease_merged[each1] = each[0]

ICD9_ICD10_map = dict()
with open('../name_map/icd10_icd9cm_map.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        ICD10 = str2[0]
        if ICD10[0] > 'N':
            continue
        if ICD10[:3] in disease_merged:
            ICD10 = disease_merged[ICD10[:3]]
        else:
            ICD10 = ICD10[:3]
        set1 = set(str2[1].split(';'))
        for each in set1:
            ICD9_ICD10_map[each] = ICD10
    infile.close()

comor_df = pd.read_csv('../phenotype/comorbidity_filter.csv', dtype=str)
all_comor_disease = set(comor_df['disease1']) | set(comor_df['disease2'])

disease_class = dict()
df = pd.read_excel('../network_random/zhouxuezhong_2018_download/mmc9.xlsx')
list1 = df.values.tolist()
for each in list1:
    category = each[0].split('.')[0]
    icd9 = each[2]
    if icd9 in ICD9_ICD10_map:
        icd10 = ICD9_ICD10_map[icd9]
        if icd10 in all_comor_disease:
            disease_class[icd10] = category
category_disease = dict()
for each in disease_class:
    category = disease_class[each]
    if category not in category_disease:
        category_disease[category] = set()
    category_disease[category].add(each)

# plot
list1 = list()
list2 = list()
category_order = list(category_disease.keys())
category_order.sort()
for each in category_order:
    list1.append(each)
    list2.append(len(category_disease[each]))

plt.bar(range(0, len(list2)), list2, color='blue')
plt.xticks(range(0, len(list1)), list1, fontsize=14, rotation=90)
plt.ylabel('Number of diseases', fontsize=14)
plt.xlabel('Category', fontsize=14)
plt.savefig('a.pdf', bbox_inches='tight')


df = pd.read_csv('../network_random/comorbid_tendency_NCD.csv')
list1 = df.values.tolist()

original = sum(list1[0][:len(category_disease)])
count = 0
for each in list1[1:]:
    s = sum(each[:len(category_disease)])
    if s >= original:
        count += 1
print('inner comorbid tendency as a whole: ' + str(float(count)/df.shape[0]))


list1 = list(df.iloc[0])
category_category_order = list(df.columns)
category_category_comorbidity_original = dict()
for i in range(0, len(list1)):
    category_category_comorbidity_original[category_category_order[i]] = list1[i]

category_category_comorbidity_significance = dict()
for each in category_category_order:
    count = df[df[each]>=df[each][0]].shape[0]
    category_category_comorbidity_significance[each] = float(count)/df.shape[0]

list1 = list()
for each in category_category_comorbidity_significance:
    c1, c2 = each.split('-')
    set1 = category_disease[c1]
    set2 = category_disease[c2]
    count = category_category_comorbidity_original[each]
    if c1 == c2:
        total = len(list(combinations(set1, 2)))
    else:
        total = len(list(product(set1, set2)))

    if total == 0:
        ratio = 0
    else:
        ratio = float(count) / total
    significance = category_category_comorbidity_significance[each]
    list1.append([c1, c2, count, ratio, significance])

df = pd.DataFrame(list1, columns=['category1', 'category2', 'comorbidity count', 'Jaccard index', 'p'])
df['q'] = multitest.fdrcorrection(df['p'])[1]
df.to_csv('../network_random/comorbidity_category_significance_NCD.txt', sep='\t', index=False)


# plot
all_category = set(df['category1']) | set(df['category2'])
category_order = list(all_category)
category_order.sort()
dict1 = dict()
i = 0
for each in category_order:
    dict1[each] = i
    i += 1

jar_data = np.zeros((len(category_order), len(category_order)))
q_data = np.zeros((len(category_order), len(category_order)))
list1 = df.values.tolist()
for each in list1:
    c1 = each[0]
    c2 = each[1]
    sig = each[5]
    jar = each[3]
    index1 = dict1[c1]
    index2 = dict1[c2]
    jar_data[index1, index2] = jar
    jar_data[index2, index1] = jar
    q_data[index1, index2] = sig
    q_data[index2, index1] = sig

df = pd.DataFrame(jar_data)
df.columns = category_order
df.index = category_order
df.to_csv('ratio.csv')

df = pd.DataFrame(q_data)
df.columns = category_order
df.index = category_order
df.to_csv('q.csv')

robjects.r.source('plot_corr.R')