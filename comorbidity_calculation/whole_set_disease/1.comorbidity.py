from scipy import stats
import numpy as np
from itertools import combinations
import pandas as pd

#--------------------------------------RR------------------------------------------#
def RR_p(Ii, Ij, Cij, N, nn):
    temp1 = Cij * N
    temp2 = Ii * Ij
    CijStar = temp2 / float(N)
    p = stats.poisson.pmf(nn, CijStar)
    p_value = 1 - sum(p[:Cij])
    if p_value < 0:
        p_value = 0
    return temp1/temp2, p_value


#--------------------------------------fai and p ------------------------------------------#
def fai_p(Ii, Ij, Cij, N):
    temp1 = N * Cij - Ii * Ij
    try:
        temp2 = np.sqrt(float(Ii * Ij * (N - Ii) * (N - Ij)))
    except:
        raise Exception
    fai = temp1/temp2
    n = max(Ii, Ij)
    temp1 = fai * np.sqrt(n-2)
    temp2 = np.sqrt(float(1-np.power(fai, 2)))
    t = temp1/temp2
    p = 1 - stats.t.cdf(t, n)
    return fai, p


if __name__ == '__main__':

    df = pd.read_table('../phenotype/merged_icd10.txt', dtype=str)
    icd10_merged = set(df['icd10'])
    dict1 = dict()
    for each in icd10_merged:
        list1 = list(each.split(';'))
        for each1 in list1:
            dict1[each1] = each

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
                if each1 not in dict1:
                    continue
                code = dict1[each1]
                if code not in disease_patient:
                    disease_patient[code] = set()
                disease_patient[code].add(person)
        infile.close()

    print(len(disease_patient))


    with open('../phenotype/comorbidity_alldisease.csv', 'w+') as outfile:
        outfile.write('disease1' + ',' + 'disease2' + ',' + 'I1' + ',' + 'I2' + ',' + 'C12' + ',' + 'RR' + ',' + 'p value' + '\n')
    outfile.close()



    N = 410309
    nn = np.arange(0, N)
    i = 0
    set1 = set(disease_patient.keys())
    all_combination = list(combinations(set1, 2))

    print('all comorbidity combination: ' + str(len(all_combination)))
    result = list()
    i = 0
    for each in all_combination:
        if i % 1000 == 0:
            print(i)
        i += 1
        code1 = each[0]
        code2 = each[1]
        patient_set1 = disease_patient[code1]
        patient_set2 = disease_patient[code2]
        comorbidity_patient = patient_set1 & patient_set2
        Ii = len(patient_set1)
        Ij = len(patient_set2)
        Cij = len(comorbidity_patient)
        rr, p = RR_p(Ii, Ij, Cij, N, nn)
        if code1 < code2:
            result.append(code1 + ',' + code2 + ',' + str(Ii) + ',' + str(Ij) + ',' + str(Cij) + ',' + str(rr) + ',' + str(p))
        else:
            result.append(code2 + ',' + code1 + ',' + str(Ij) + ',' + str(Ii) + ',' + str(Cij) + ',' + str(rr) + ',' + str(p))
        if len(result) == 3000:
            with open('../phenotype/comorbidity_alldisease.csv', 'a') as outfile:
                outfile.write('\n'.join(result) + '\n')
            outfile.close()
            result = list()
    with open('../phenotype/comorbidity_alldisease.csv', 'a') as outfile:
        outfile.write('\n'.join(result) + '\n')
    outfile.close()