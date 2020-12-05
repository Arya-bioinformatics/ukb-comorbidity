import pandas as pd
from itertools import combinations
from statsmodels.stats import multitest
from itertools import product
import numpy as np
import rpy2.robjects as robjects


category_disease = dict()
disease_category = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    code = each[0]
    category = each[2]
    disease_category[code] = category
    if category not in category_disease:
        category_disease[category] = set()
    category_disease[category].add(code)


df = pd.read_csv('../network_random/comorbid_tendency.csv')
list1 = df.values.tolist()

original = sum(list1[0][:24])
count = 0
for each in list1[1:]:
    s = sum(each[:24])
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
df.to_csv('../network_random/comorbidity_category_significance.txt', sep='\t', index=False)


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