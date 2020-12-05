import pandas as pd
from itertools import combinations
from datetime import datetime


icd10_merged = dict()
disease_category = dict()
disease_code_description = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    for icd10 in each[0].split(';'):
        icd10_merged[icd10] = each[0]
    disease_category[each[0]] = each[2]
    disease_code_description[each[0]] = each[1]


patient_disease = dict()
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
                each1 = icd10_merged[each[:3]]
            else:
                each1 = each[:3]
            list2.append(each1)
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
        patient_diagnosis_date[person] = list1
    infile.close()

# ----------------------- same date --------------------------- #
i = 0
double_dis_pat = dict()
for each in patient_disease:
    i += 1
    if i % 1000 == 0:
        print(i)
    diagnosis = patient_disease[each]
    if len(set(diagnosis)) < 2:
        continue
    date = patient_diagnosis_date[each]
    if len(date) == len(set(date)):
        continue
    date_count = pd.value_counts(date)
    duplicated_date = set(date_count[date_count>1].index)
    df = pd.DataFrame([diagnosis, date])
    df = df.T
    df.columns = ['diagnosis', 'date']
    df = df.drop_duplicates(['diagnosis', 'date'])
    df = df[df['date'].isin(duplicated_date)]
    all_date = set(df['date'])
    for each_date in all_date:
        these_diseases = df.loc[df['date']==each_date, 'diagnosis']
        set1 = set(these_diseases)
        set2 = set()
        for each1 in set1:
            if each1[0] <= 'N':
                set2.add(each1)
        if len(set2) < 2:
            continue
        all_combination = combinations(set2, 2)
        for each2 in all_combination:
            code1 = each2[0]
            code2 = each2[1]
            if code1 < code2:
                double_disease = code1 + '-' + code2
            else:
                double_disease = code2 + '-' + code1
            if double_disease not in double_dis_pat:
                double_dis_pat[double_disease] = set()
            double_dis_pat[double_disease].add(each)
with open('../phenotype/double_disease_patient_phecode.txt', 'w+') as outfile:
    outfile.write('comorbidity_candidate' + '\t' + 'patient' + '\n')
    for each in double_dis_pat:
        outfile.write(each + '\t' + ';'.join(double_dis_pat[each]) + '\n')
    outfile.close()

# ----------------------- half year --------------------------- #
i = 0
double_dis_pat = dict()
for each in patient_disease:
    i += 1
    if i % 1000 == 0:
        print(i)
    diagnosis = patient_disease[each]
    if len(set(diagnosis)) < 2:
        continue
    date = patient_diagnosis_date[each]
    df = pd.DataFrame([diagnosis, date])
    df = df.T
    df.columns = ['diagnosis', 'date']
    for each1 in combinations(date, 2):
        y1, m1, d1 = each1[0].split('-')
        y2, m2, d2 = each1[1].split('-')
        delta = datetime(int(y1), int(m1), int(d1)) - datetime(int(y2), int(m2), int(d2))
        if abs(delta.days) < 365/2:
            list3 = list(df.loc[df['date'].isin(set(each1)), 'diagnosis'])
            for code1, code2 in combinations(list3, 2):
                if (code1[0] > 'N') | (code2[0] > 'N'):
                    continue
                if code1 == code2:
                    continue
                if code1 < code2:
                    key = code1 + '-' + code2
                else:
                    key = code2 + '-' + code1
                if key not in double_dis_pat:
                    double_dis_pat[key] = set()
                double_dis_pat[key].add(each)


with open('../phenotype/double_disease_patient_phecode_0.5year.txt', 'w+') as outfile:
    outfile.write('comorbidity_candidate' + '\t' + 'patient' + '\n')
    for each in double_dis_pat:
        outfile.write(each + '\t' + ';'.join(double_dis_pat[each]) + '\n')
    outfile.close()


# ----------------------- one year --------------------------- #
i = 0
double_dis_pat = dict()
for each in patient_disease:
    i += 1
    if i % 1000 == 0:
        print(i)
    diagnosis = patient_disease[each]
    if len(set(diagnosis)) < 2:
        continue
    date = patient_diagnosis_date[each]
    df = pd.DataFrame([diagnosis, date])
    df = df.T
    df.columns = ['diagnosis', 'date']
    for each1 in combinations(date, 2):
        y1, m1, d1 = each1[0].split('-')
        y2, m2, d2 = each1[1].split('-')
        delta = datetime(int(y1), int(m1), int(d1)) - datetime(int(y2), int(m2), int(d2))
        if abs(delta.days) < 365:
            list3 = list(df.loc[df['date'].isin(set(each1)), 'diagnosis'])
            for code1, code2 in combinations(list3, 2):
                if (code1[0] > 'N') | (code2[0] > 'N'):
                    continue
                if code1 == code2:
                    continue
                if code1 < code2:
                    key = code1 + '-' + code2
                else:
                    key = code2 + '-' + code1
                if key not in double_dis_pat:
                    double_dis_pat[key] = set()
                double_dis_pat[key].add(each)


with open('../phenotype/double_disease_patient_phecode_1year.txt', 'w+') as outfile:
    outfile.write('comorbidity_candidate' + '\t' + 'patient' + '\n')
    for each in double_dis_pat:
        outfile.write(each + '\t' + ';'.join(double_dis_pat[each]) + '\n')
    outfile.close()


# ----------------------- 2 years --------------------------- #
i = 0
double_dis_pat = dict()
for each in patient_disease:
    i += 1
    # if i % 1000 == 0:
    #     print(i)
    print(i)
    diagnosis = patient_disease[each]
    if len(set(diagnosis)) < 2:
        continue
    date = patient_diagnosis_date[each]
    df = pd.DataFrame([diagnosis, date])
    df = df.T
    df.columns = ['diagnosis', 'date']
    for each1 in combinations(date, 2):
        y1, m1, d1 = each1[0].split('-')
        y2, m2, d2 = each1[1].split('-')
        delta = datetime(int(y1), int(m1), int(d1)) - datetime(int(y2), int(m2), int(d2))
        if abs(delta.days) < 365*2:
            list3 = list(df.loc[df['date'].isin(set(each1)), 'diagnosis'])
            for code1, code2 in combinations(list3, 2):
                if (code1[0] > 'N') | (code2[0] > 'N'):
                    continue
                if code1 == code2:
                    continue
                if code1 < code2:
                    key = code1 + '-' + code2
                else:
                    key = code2 + '-' + code1
                if key not in double_dis_pat:
                    double_dis_pat[key] = set()
                double_dis_pat[key].add(each)

with open('../phenotype/double_disease_patient_phecode_2year.txt', 'w+') as outfile:
    outfile.write('comorbidity_candidate' + '\t' + 'patient' + '\n')
    for each in double_dis_pat:
        outfile.write(each + '\t' + ';'.join(double_dis_pat[each]) + '\n')
    outfile.close()

# ----------------------- 3 years --------------------------- #
i = 0
double_dis_pat = dict()
for each in patient_disease:
    i += 1
    if i % 1000 == 0:
        print(i)
    diagnosis = patient_disease[each]
    if len(set(diagnosis)) < 2:
        continue
    date = patient_diagnosis_date[each]
    df = pd.DataFrame([diagnosis, date])
    df = df.T
    df.columns = ['diagnosis', 'date']
    for each1 in combinations(date, 2):
        y1, m1, d1 = each1[0].split('-')
        y2, m2, d2 = each1[1].split('-')
        delta = datetime(int(y1), int(m1), int(d1)) - datetime(int(y2), int(m2), int(d2))
        if abs(delta.days) < 365*3:
            list3 = list(df.loc[df['date'].isin(set(each1)), 'diagnosis'])
            for code1, code2 in combinations(list3, 2):
                if (code1[0] > 'N') | (code2[0] > 'N'):
                    continue
                if code1 == code2:
                    continue
                if code1 < code2:
                    key = code1 + '-' + code2
                else:
                    key = code2 + '-' + code1
                if key not in double_dis_pat:
                    double_dis_pat[key] = set()
                double_dis_pat[key].add(each)


with open('../phenotype/double_disease_patient_phecode_3year.txt', 'w+') as outfile:
    outfile.write('comorbidity_candidate' + '\t' + 'patient' + '\n')
    for each in double_dis_pat:
        outfile.write(each + '\t' + ';'.join(double_dis_pat[each]) + '\n')
    outfile.close()


# ----------------------- 4 years --------------------------- #
i = 0
double_dis_pat = dict()
for each in patient_disease:
    i += 1
    if i % 1000 == 0:
        print(i)
    diagnosis = patient_disease[each]
    if len(set(diagnosis)) < 2:
        continue
    date = patient_diagnosis_date[each]
    df = pd.DataFrame([diagnosis, date])
    df = df.T
    df.columns = ['diagnosis', 'date']
    for each1 in combinations(date, 2):
        y1, m1, d1 = each1[0].split('-')
        y2, m2, d2 = each1[1].split('-')
        delta = datetime(int(y1), int(m1), int(d1)) - datetime(int(y2), int(m2), int(d2))
        if abs(delta.days) < 365*4:
            list3 = list(df.loc[df['date'].isin(set(each1)), 'diagnosis'])
            for code1, code2 in combinations(list3, 2):
                if (code1[0] > 'N') | (code2[0] > 'N'):
                    continue
                if code1 == code2:
                    continue
                if code1 < code2:
                    key = code1 + '-' + code2
                else:
                    key = code2 + '-' + code1
                if key not in double_dis_pat:
                    double_dis_pat[key] = set()
                double_dis_pat[key].add(each)


with open('../phenotype/double_disease_patient_phecode_4year.txt', 'w+') as outfile:
    outfile.write('comorbidity_candidate' + '\t' + 'patient' + '\n')
    for each in double_dis_pat:
        outfile.write(each + '\t' + ';'.join(double_dis_pat[each]) + '\n')
    outfile.close()


# ----------------------- 5 years --------------------------- #
i = 0
double_dis_pat = dict()
for each in patient_disease:
    i += 1
    if i % 1000 == 0:
        print(i)
    diagnosis = patient_disease[each]
    if len(set(diagnosis)) < 2:
        continue
    date = patient_diagnosis_date[each]
    df = pd.DataFrame([diagnosis, date])
    df = df.T
    df.columns = ['diagnosis', 'date']
    for each1 in combinations(date, 2):
        y1, m1, d1 = each1[0].split('-')
        y2, m2, d2 = each1[1].split('-')
        delta = datetime(int(y1), int(m1), int(d1)) - datetime(int(y2), int(m2), int(d2))
        if abs(delta.days) < 365*5:
            list3 = list(df.loc[df['date'].isin(set(each1)), 'diagnosis'])
            for code1, code2 in combinations(list3, 2):
                if (code1[0] > 'N') | (code2[0] > 'N'):
                    continue
                if code1 == code2:
                    continue
                if code1 < code2:
                    key = code1 + '-' + code2
                else:
                    key = code2 + '-' + code1
                if key not in double_dis_pat:
                    double_dis_pat[key] = set()
                double_dis_pat[key].add(each)

with open('../phenotype/double_disease_patient_phecode_5year.txt', 'w+') as outfile:
    outfile.write('comorbidity_candidate' + '\t' + 'patient' + '\n')
    for each in double_dis_pat:
        outfile.write(each + '\t' + ';'.join(double_dis_pat[each]) + '\n')
    outfile.close()