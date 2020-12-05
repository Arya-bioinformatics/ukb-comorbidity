import pandas as pd
from scipy import stats


icd10_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    code = each[0]
    list2 = list(code.split(';'))
    for each1 in list2:
        icd10_merged[each1] = code


patient_disease = dict()
disease_patient = dict()
with open('../phenotype/field41270.csv', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str1 = str1.replace('"', '')
        str2 = str1.split(',')
        person = str2[0]
        list1 = str2[1:]
        while '' in list1:
            list1.remove('')
        if len(list1) == 0:
            continue
        list2 = list()
        for each in list1:
            if each[:3] in icd10_merged:
                list2.append(icd10_merged[each[:3]])
                if icd10_merged[each[:3]] not in disease_patient:
                    disease_patient[icd10_merged[each[:3]]] = set()
                disease_patient[icd10_merged[each[:3]]].add(person)
            else:
                list2.append(each)
                if each not in disease_patient:
                    disease_patient[each] = set()
                disease_patient[each].add(person)
        patient_disease[person] = list2
    infile.close()


patient_diagnosis_date = dict()
with open('../phenotype/field41280.csv', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str1 = str1.replace('"', '')
        str2 = str1.split(',')
        person = str2[0]
        list1 = str2[1:]
        while '' in list1:
            list1.remove('')
        if len(list1) == 0:
            continue
        list2 = list()
        for each in list1:
            list2.append(int(''.join(each.split('-'))))
        patient_diagnosis_date[person] = list2
    infile.close()

double_disease_patient = dict()
with open('../phenotype/double_disease_patient_phecode.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient[str2[0]] = set(str2[1].split(';'))
    infile.close()


df = pd.read_csv('../phenotype/comorbidity_filter.csv')
list1 = df.values.tolist()
result = list()
direct_dict = dict()
i = 0
for each in list1:
    if i % 100 == 0:
        print(i)
    i += 1
    code1, code2 = each[0], each[1]
    comor_person = disease_patient[code1] & disease_patient[code2]
    Na = 0
    Nb = 0
    double_patients = double_disease_patient[code1 + '-' + code2]
    Nsame = len(double_patients)
    for each1 in comor_person:
        if each1 in double_patients:
            continue
        df1 = pd.DataFrame()
        df1['disease'] = patient_disease[each1]
        df1['date'] = patient_diagnosis_date[each1]
        date1 = min(df1.loc[df1['disease']==code1, 'date'])
        date2 = min(df1.loc[df1['disease']==code2, 'date'])
        if date1 < date2:
            Na += 1
        elif date1 > date2:
            Nb += 1
        else:
            raise Exception
    pa = stats.binom_test(Na, Na + Nb + Nsame, p=0.5, alternative='greater')
    pb = stats.binom_test(Nb, Na + Nb + Nsame, p=0.5, alternative='greater')
    result.append([code1, code2, pa])
    result.append([code2, code1, pb])
    direct_dict[code1 + '-' + code2] = [Na, Nb, pa, pb]
df1 = pd.DataFrame(result, columns=['from', 'to', 'p'])
df1 = df1[df1['p']<0.05/df1.shape[0]]
list2 = df1.values.tolist()
set1 = set()
for each in list2:
    set1.add(each[0] + '-' + each[1])
result = list()
for each in list1:
    if each[0] + '-' + each[1] in set1:
        direction = 1
    elif each[1] + '-' + each[0] in set1:
        direction = -1
    else:
        direction = 0
    if each[0] + '-' + each[1] not in direct_dict:
        continue
    result.append(each + [direction] + direct_dict[each[0] + '-' + each[1]])
df = pd.DataFrame(result, columns=list(df.columns) + ['direction', 'N1', 'N2', 'P1_direction', 'p2_direction'])
df.to_csv('../phenotype/comorbidity_filter_direction.csv', index=False)