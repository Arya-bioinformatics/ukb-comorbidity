import pandas as pd
from scipy.stats import pearsonr

disease_gene = dict()
with open('../genome/disease_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease_gene[str2[0]] = len(str2[1].split(';'))
    infile.close()


disease_comorbid_gene = dict()
with open('../overlap/comorbidity_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1 = str2[0]
        code2 = str2[1]
        set1 = set(str2[4].split(';'))
        if code1 not in disease_comorbid_gene:
            disease_comorbid_gene[code1] = set()
        disease_comorbid_gene[code1] |= set1
        if code2 not in disease_comorbid_gene:
            disease_comorbid_gene[code2] = set()
        disease_comorbid_gene[code2] |= set1
    infile.close()


disease_comorbidity = dict()
df = pd.read_csv('../phenotype/comorbidity_filter.csv')
all_disease = set(df['disease1']) | set(df['disease2'])
for each in all_disease:
    disease_comorbidity[each] = df[(df['disease1'] == each) | (df['disease2'] == each)].shape[0]

list1 = list()
list2 = list()
for each in disease_comorbidity:
    if each in disease_gene:
        list1.append(disease_comorbidity[each])
        list2.append(disease_gene[each])

[cof, p] = pearsonr(list1, list2)
print([cof, p])


list1 = list()
list2 = list()
for each in disease_comorbidity:
    if each in disease_comorbid_gene:
        list1.append(disease_comorbidity[each])
        list2.append(len(disease_comorbid_gene[each]))

[cof, p] = pearsonr(list1, list2)
print([cof, p])