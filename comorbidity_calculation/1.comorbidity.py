import pandas as pd
from scipy import stats
import numpy as np
from itertools import combinations



def RR(Ii, Ij, Cij, N):
    temp1 = Cij * N
    temp2 = Ii * Ij
    return float(temp1)/float(temp2)

def significance(Ii, Ij, Cij, N, nn):
    CijStar = float(Ii * Ij)/float(N)
    p = stats.poisson.pmf(nn, CijStar)
    p_value = 1 - sum(p[:Cij])
    if p_value < 0:
        p_value = 0
    return p_value


icd10_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    for icd10 in each[0].split(';'):
        icd10_merged[icd10] = each[0]

# disease - patient
disease_patient = dict()
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
        for each in set1:
            each1 = each[:3]
            if each1 in icd10_merged:
                disease = icd10_merged[each1]
                if disease not in disease_patient:
                    disease_patient[disease] = set()
                disease_patient[disease].add(person)
    infile.close()


with open('../phenotype/comorbidity.csv', 'w+') as comorbidityFile:
    comorbidityFile.write('disease1' + ',' + 'disease2' + ',' + 'patient number of disease1' + ',' + 'patient number of disease2' + ',' + 'patiant number of disease1 and disease2' + ','
                              + 'RR' + ',' + 'p_value' + '\n')
    comorbidityFile.close()


N = 410293
nn = np.arange(0, N)
i = 0
set1 = set(disease_patient.keys())
all_combination = list(combinations(set1, 2))

print('all comorbidity combination: ' + str(len(all_combination)))
result = list()
for each in all_combination:
    if i % 1000 == 0:
        print(i)
    i = i + 1
    code1 = each[0]
    code2 = each[1]
    patient_set1 = disease_patient[code1]
    patient_set2 = disease_patient[code2]
    comorbidity_patient = patient_set1 & patient_set2
    Ii = len(patient_set1)
    Ij = len(patient_set2)
    Cij = len(comorbidity_patient)
    rr = RR(Ii, Ij, Cij, N)
    p_value = significance(Ii, Ij, Cij, N, nn)
    if code1 < code2:
        result.append(code1 + ',' + code2 + ',' + str(Ii) + ',' + str(Ij) + ',' + str(Cij) + ',' + str(rr) + ',' + str(p_value))
    else:
        result.append(code2 + ',' + code1 + ',' + str(Ij) + ',' + str(Ii) + ',' + str(Cij) + ',' + str(rr) + ',' + str(p_value))
    if len(result) == 1000:
        with open('../phenotype/comorbidity.csv', 'a') as comorbidityFile:
            comorbidityFile.write('\n'.join(result) + '\n')
        comorbidityFile.close()
        result = list()

with open('../phenotype/comorbidity.csv', 'a') as comorbidityFile:
    comorbidityFile.write('\n'.join(result) + '\n')
    comorbidityFile.close()