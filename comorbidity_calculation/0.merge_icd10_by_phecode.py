import pandas as pd


disease_code_description = dict()
with open('../phenotype/coding19.tsv', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code = str2[0]
        des = str2[1]
        if len(code) == 3:
            disease_code_description[code] = des
    infile.close()

all_code = set()
with open('../phenotype/field41270.csv', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str1 = str1.replace('"', '')
        str2 = str1.split(',')
        set1 = set(str2[1:])
        if '' in set1:
            set1.remove('')
        if len(set1) == 0:
            continue
        for each in set1:
            each1 = each[:3]
            if each1[0] > 'N':
                continue
            all_code.add(each1)
    infile.close()

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
            if each1[0] > 'N':
                continue
            if each1 not in disease_patient:
                disease_patient[each1] = set()
            disease_patient[each1].add(person)
    infile.close()

common_disease = set()
for each in disease_patient:
    if len(disease_patient[each]) < 410:
        continue
    if (each == 'C77') | (each == 'C78') | (each == 'C79'):
        continue
    common_disease.add(each)

df = pd.read_csv('../phenotype/phecode_icd10.csv', dtype=str)
df = df.dropna(axis=0, how='any')
list1 = df[['ICD10', 'PheCode']].values.tolist()
ICD10_phecode = dict()
for each in list1:
    icd10 = each[0]
    phecode = each[1]
    if '.' in phecode:
        [s1, s2] = phecode.split('.')
        phecode = s1 + '.' + s2[0]
    ICD10_phecode[icd10] = phecode

list1 = list()
for each in common_disease:
    if each in ICD10_phecode:
        phecode = ICD10_phecode[each]
    else:
        phecode = each
    list1.append([each, phecode])

df = pd.DataFrame(list1, columns=['icd10', 'phecode'])
all_phecode = set(df['phecode'])
result = list()
for each in all_phecode:
    list1 = list(df.loc[df['phecode']==each, 'icd10'])
    list1.sort()
    list2 = list()
    for each1 in list1:
        des = disease_code_description[each1]
        list2.append(des)
    result.append([';'.join(list1), ';'.join(list2)])

df = pd.DataFrame(result, columns=['ICD10', 'Description'])
df.to_csv('../phenotype/merged_icd10.csv', index=False)