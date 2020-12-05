import pandas as pd
import networkx as nx
from itertools import combinations
import random


disease_class = dict()
disease_index = dict()
disease_index1 = dict()
all_category = set()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
i = 0
for each in list1:
    code = each[0]
    category = each[2]
    disease_class[code] = category
    disease_index[i] = code
    disease_index1[code] = i
    i += 1
    all_category.add(category)
all_category = list(all_category)

set1 = set()
df = pd.read_csv('../phenotype/comorbidity_filter.csv')
list1 = df.values.tolist()
for each in list1:
    code1 = each[0]
    code2 = each[1]
    index1 = disease_index1[code1]
    index2 = disease_index1[code2]
    set1.add((index1, index2))

G = nx.Graph()
G.add_edges_from(set1)
list1 = list(G.edges)

all_category.sort()
category_category_order = list()
for each in all_category:
    category_category_order.append(each + '-' + each)
comb = list(combinations(all_category, 2))
for each in comb:
    (c1, c2) = each
    if c1 < c2:
        category_category_order.append((c1 + '-' + c2))
    else:
        category_category_order.append((c2 + '-' + c1))

result = list()
i = 0
while i < 1000001:
    print(i)
    this_result = dict()
    for each in list1:
        [index1, index2] = each
        code1 = disease_index[index1]
        code2 = disease_index[index2]
        c1 = disease_class[code1]
        c2 = disease_class[code2]
        if c1 < c2:
            key = c1 + '-' + c2
        else:
            key = c2 + '-' + c1
        if key not in this_result:
            this_result[key] = 0
        this_result[key] += 1
    temp = list()
    for each in category_category_order:
        if each in this_result:
            count = this_result[each]
        else:
            count = 0
        temp.append(count)
    result.append(temp)

    temp1 = list(disease_index.keys())
    temp2 = list(disease_index.values())
    j = 0
    while j < 10:
        random.shuffle(temp1)
        random.shuffle(temp2)
        j += 1
    disease_index = dict()
    for m in range(0, len(temp1)):
        disease_index[temp1[m]] = temp2[m]
    i += 1

df = pd.DataFrame(result, columns=category_category_order)
df.to_csv('../network_random/comorbid_tendency.csv', index=False)