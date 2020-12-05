import pandas as pd
from itertools import combinations, product
from scipy.stats import fisher_exact
import numpy as np


disease_merged = dict()
disease_description = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    set1 = set(each[0].split(';'))
    disease_description[each[0]] = each[1]
    for each1 in set1:
        disease_merged[each1] = each[0]

# same day
double_disease_patient = dict()
with open('../phenotype/double_disease_patient_phecode.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient[str2[0]] = len(set(str2[1].split(';')))
    infile.close()

disease_pairs_prefiltered = set()
df = pd.read_csv('../phenotype/comorbidity.csv')
list1 = df.values.tolist()
for each in list1:
    code1, code2, Ii, Ij = each[:4]
    if code1 + '-' + code2 in double_disease_patient:
        count = double_disease_patient[code1 + '-' + code2]
    if (count/Ii>0.01) | (count/Ij>0.01):
        disease_pairs_prefiltered.add((code1, code2))


# ----------------- compare ukb and hidalgo ----------------- #
ICD9_ICD10_map = dict()
with open('../name_map/icd10_icd9cm_map.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        ICD10 = str2[0]
        if len(ICD10) != 3:
            continue
        if ICD10[0] > 'N':
            continue
        if ICD10 in disease_merged:
            ICD10 = disease_merged[ICD10]
        set1 = set(str2[1].split(';'))
        for each in set1:
            if len(each) != 3:
                continue
            ICD9_ICD10_map[each] = ICD10
    infile.close()

df1 = pd.read_csv('../phenotype/comorbidity_filter.csv')
ukb_disease = set(df1['disease1']) | set(df1['disease2'])


df = pd.read_table('../phenotype/comorbidity_Hidalgo/AllNet3.tsv')
# bonferroni correction
df_bf = df[df['RR']>1]
df_bf = df_bf[df_bf['pval']<0.05/df_bf.shape[0]]
list1 = df_bf.values.tolist()
for each in list1:
    d1, d2 = each[0], each[1]
    temp1 = ('00' + str(int(d1)))[-3:]
    temp2 = ('00' + str(int(d2)))[-3:]
    if (temp1 in ICD9_ICD10_map) & (temp2 in ICD9_ICD10_map):
        code1 = ICD9_ICD10_map[temp1]
        code2 = ICD9_ICD10_map[temp2]
    else:
        continue
    if (code1 in ukb_disease) & (code2 in ukb_disease):
        print(disease_description[code1] + '-' + disease_description[code2]) # only one with two diseases in ukb diseases

df2 = df[(df['RR']>1) & (df['pval.adj']<0.05)]
hidalgo_disease = set(df2['#icd9-1']) | set(df2['icd9-2'])


ukb_hidalgo_common_disease = set()
for each in hidalgo_disease:
    temp = '00' + str(each)
    temp = temp[-3:]
    if temp in ICD9_ICD10_map:
        ICD10 = ICD9_ICD10_map[temp]
        if ICD10 in ukb_disease:
            ukb_hidalgo_common_disease.add(ICD10)


ukb_comorbidity = set()
list1 = df1.values.tolist()
for each in list1:
    code1, code2 = each[0], each[1]
    if code1 not in ukb_hidalgo_common_disease:
        continue
    if code2 not in ukb_hidalgo_common_disease:
        continue
    ukb_comorbidity.add(code1 + '-' + code2)

hidalgo_comorbidity = set()
list2 = df2.values.tolist()
for each in list2:
    temp1, temp2 = int(each[0]), int(each[1])
    temp1 = ('00' + str(temp1))[-3:]
    temp2 = ('00' + str(temp2))[-3:]
    if temp1 not in ICD9_ICD10_map:
        continue
    if temp2 not in ICD9_ICD10_map:
        continue
    code1 = ICD9_ICD10_map[temp1]
    code2 = ICD9_ICD10_map[temp2]
    if (code1 in ukb_hidalgo_common_disease) & (code2 in ukb_hidalgo_common_disease):
        if code1 < code2:
            hidalgo_comorbidity.add(code1 + '-' + code2)
        else:
            hidalgo_comorbidity.add(code2 + '-' + code1)

print('------------------ compare ukb and hidalgo -----------------')
print('ukb and hidalgo common diseases: ' + str(len(ukb_hidalgo_common_disease)))
print('involved ukb comorbidity: ' + str(len(ukb_comorbidity)))
print('involved hidalgo comorbidity: ' + str(len(hidalgo_comorbidity)))
print('intersaction: ' + str(len(ukb_comorbidity & hidalgo_comorbidity)))
print(ukb_comorbidity & hidalgo_comorbidity)

a = len(ukb_comorbidity & hidalgo_comorbidity)
b = len(hidalgo_comorbidity) - a
c = len(ukb_comorbidity) - a
d = len(set(combinations(ukb_hidalgo_common_disease, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])


# ----------------- compare ukb and jensen ------------------ #
jensen_disease = set()
df3 = pd.read_excel('../phenotype/comorbidity_Jensen/ncomms5022-s2.xlsx')
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    jensen_disease.add(code1)
    jensen_disease.add(code2)

ukb_jensen_common_disease = ukb_disease & jensen_disease

ukb_comorbidity = set()
list1 = df1.values.tolist()
for each in list1:
    code1, code2 = each[0], each[1]
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    ukb_comorbidity.add(code1 + '-' + code2)

jensen_comorbidity = set()
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    if code1 < code2:
        jensen_comorbidity.add(code1 + '-' + code2)
    else:
        jensen_comorbidity.add(code2 + '-' + code1)

print('------------------ compare ukb and jensen -----------------')
print('ukb and jensen common diseases: ' + str(len(ukb_jensen_common_disease)))
print('involved ukb comorbidity: ' + str(len(ukb_comorbidity)))
print('involved jensen comorbidity: ' + str(len(jensen_comorbidity)))
print('intersaction: ' + str(len(ukb_comorbidity & jensen_comorbidity)))

a = len(ukb_comorbidity & jensen_comorbidity)
b = len(jensen_comorbidity) - a
c = len(ukb_comorbidity) - a
d = len(set(combinations(ukb_jensen_common_disease, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])

# ----------------- compare directional ukb and jensen ------------------ #
jensen_disease = set()
df3 = pd.read_excel('../phenotype/comorbidity_Jensen/ncomms5022-s2.xlsx')
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    jensen_disease.add(code1)
    jensen_disease.add(code2)

df4 = pd.read_csv('../phenotype/comorbidity_filter_direction.csv')
ukb_directional_disease = set(df4.loc[df4['direction']!=0, 'disease1']) | set(df4.loc[df4['direction']!=0, 'disease2'])

ukb_jensen_common_disease = ukb_directional_disease & jensen_disease

ukb_comorbidity = set()
list1 = df4.values.tolist()
for each in list1:
    code1, code2, direction = each[0], each[1], each[12]
    if direction == 0:
        continue
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    ukb_comorbidity.add(code1 + '-' + code2)

jensen_comorbidity = set()
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    if code1 < code2:
        jensen_comorbidity.add(code1 + '-' + code2)
    else:
        jensen_comorbidity.add(code2 + '-' + code1)

print('------------------ compare directional ukb and jensen -----------------')
print('ukb and jensen common diseases: ' + str(len(ukb_jensen_common_disease)))
print('involved ukb comorbidity: ' + str(len(ukb_comorbidity)))
print('involved jensen comorbidity: ' + str(len(jensen_comorbidity)))
print('intersaction: ' + str(len(ukb_comorbidity & jensen_comorbidity)))

a = len(ukb_comorbidity & jensen_comorbidity)
b = len(jensen_comorbidity) - a
c = len(ukb_comorbidity) - a
d = len(set(combinations(ukb_jensen_common_disease, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])



# ----------------- compare directional-same day ukb and hidalgo ------------------ #
jensen_disease = set()
df3 = pd.read_excel('../phenotype/comorbidity_Jensen/ncomms5022-s2.xlsx')
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    jensen_disease.add(code1)
    jensen_disease.add(code2)

df4 = pd.read_csv('../phenotype/comorbidity_filter_direction.csv')
ukb_directional_disease = set(df4.loc[df4['direction']!=0, 'disease1']) | set(df4.loc[df4['direction']!=0, 'disease2'])

ukb_jensen_common_disease = ukb_directional_disease & jensen_disease

ukb_comorbidity = set()
list1 = df4.values.tolist()
for each in list1:
    code1, code2, direction = each[0], each[1], each[12]
    if direction == 0:
        continue
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    ukb_comorbidity.add(code1 + '-' + code2)

jensen_comorbidity = set()
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    if code1 < code2:
        if (code1, code2) in disease_pairs_prefiltered:
            jensen_comorbidity.add(code1 + '-' + code2)
    else:
        if (code2, code1) in disease_pairs_prefiltered:
            jensen_comorbidity.add(code2 + '-' + code1)

print('------------------ compare directional-same day ukb and jensen -----------------')
print('ukb and jensen common diseases: ' + str(len(ukb_jensen_common_disease)))
print('involved ukb comorbidity: ' + str(len(ukb_comorbidity)))
print('involved jensen comorbidity: ' + str(len(jensen_comorbidity)))
print('intersaction: ' + str(len(ukb_comorbidity & jensen_comorbidity)))

a = len(ukb_comorbidity & jensen_comorbidity)
b = len(jensen_comorbidity) - a
c = len(ukb_comorbidity) - a
set1 = set()
for code1, code2 in combinations(ukb_jensen_common_disease, 2):
    if code1 < code2:
        if (code1, code2) in disease_pairs_prefiltered:
            set1.add((code1, code2))
    else:
        if (code2, code1) in disease_pairs_prefiltered:
            set1.add((code2, code1))
d = len(set1) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])

diagnosis_count = list()
with open('../phenotype/field41270.csv', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str1 = str1.replace('"', '')
        str2 = str1.split(',')
        person = str2[0]
        set1 = set(str2[1:])
        if '' in set1:
            set1.remove('')
        if len(set1) == 0:
            continue
        set2 = set()
        for each in set1:
            set2.add(each[:3])
        diagnosis_count.append(len(set2))
    infile.close()

print('average level 2 diagnosis per person: ' + str(np.average(diagnosis_count)))