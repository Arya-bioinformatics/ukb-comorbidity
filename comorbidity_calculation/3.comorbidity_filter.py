import pandas as pd
from matplotlib import pyplot as plt


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

# 1 year
double_disease_patient_1year = dict()
with open('../phenotype/double_disease_patient_phecode_1year.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient_1year[str2[0]] = len(set(str2[1].split(';')))
    infile.close()

# 2 year
double_disease_patient_2year = dict()
with open('../phenotype/double_disease_patient_phecode_2year.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient_2year[str2[0]] = len(set(str2[1].split(';')))
    infile.close()

# 3 year
double_disease_patient_3year = dict()
with open('../phenotype/double_disease_patient_phecode_3year.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient_3year[str2[0]] = len(set(str2[1].split(';')))
    infile.close()

# 4 year
double_disease_patient_4year = dict()
with open('../phenotype/double_disease_patient_phecode_4year.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient_4year[str2[0]] = len(set(str2[1].split(';')))
    infile.close()

# 5 year
double_disease_patient_5year = dict()
with open('../phenotype/double_disease_patient_phecode_5year.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient_5year[str2[0]] = len(set(str2[1].split(';')))
    infile.close()

# half year year
double_disease_patient_halfyear = dict()
with open('../phenotype/double_disease_patient_phecode_0.5year.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        double_disease_patient_halfyear[str2[0]] = len(set(str2[1].split(';')))
    infile.close()


df = pd.read_csv('../phenotype/comorbidity.csv')
list1 = df.values.tolist()
list2 = list()
list3 = list()
list4 = list()
list5 = list()
list6 = list()
list7 = list()
list8 = list()
for each in list1:
    code1 = each[0]
    code2 = each[1]
    category1 = disease_category[code1]
    category2 = disease_category[code2]
    c1 = int(each[2])
    c2 = int(each[3])
    c12 = int(each[4])
    rr = float(each[5])
    p = float(each[6])
    name1 = disease_code_description[code1]
    name2 = disease_code_description[code2]
    key = code1 + '-' + code2
    if key not in double_disease_patient:
        continue
    cday = double_disease_patient[key]
    chalfyear = double_disease_patient_halfyear[key]
    c1year = double_disease_patient_1year[key]
    c2year = double_disease_patient_2year[key]
    c3year = double_disease_patient_3year[key]
    c4year = double_disease_patient_4year[key]
    c5year = double_disease_patient_5year[key]
    list2.append([code1, code2, name1, name2, c1, c2, c12, cday, float(cday) / c1, float(cday) / c2, rr, p, category1, category2])
    list3.append([code1, code2, name1, name2, c1, c2, c12, c1year, float(c1year) / c1, float(c1year) / c2, rr, p, category1, category2])
    list4.append([code1, code2, name1, name2, c1, c2, c12, c2year, float(c2year) / c1, float(c2year) / c2, rr, p, category1, category2])
    list5.append([code1, code2, name1, name2, c1, c2, c12, c3year, float(c3year) / c1, float(c3year) / c2, rr, p, category1, category2])
    list6.append([code1, code2, name1, name2, c1, c2, c12, c4year, float(c4year) / c1, float(c4year) / c2, rr, p, category1, category2])
    list7.append([code1, code2, name1, name2, c1, c2, c12, c5year, float(c5year) / c1, float(c5year) / c2, rr, p, category1, category2])
    list8.append([code1, code2, name1, name2, c1, c2, c12, chalfyear, float(chalfyear) / c1, float(chalfyear) / c2, rr, p, category1, category2])


list1 = list()
df1 = pd.DataFrame(list2, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[(df1['RR'] > 1) & ((df1['ratio1'] >= 0.01) | (df1['ratio2'] >= 0.01))]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
df1.to_csv('../phenotype/comorbidity_filter.csv', index=False)
list1.append(df1.shape[0])

df1 = pd.DataFrame(list8, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[(df1['RR'] > 1) & ((df1['ratio1'] >= 0.01) | (df1['ratio2'] >= 0.01))]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
list1.append(df1.shape[0])

df1 = pd.DataFrame(list3, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[(df1['RR'] > 1) & ((df1['ratio1'] >= 0.01) | (df1['ratio2'] >= 0.01))]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
list1.append(df1.shape[0])

df1 = pd.DataFrame(list4, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[(df1['RR'] > 1) & ((df1['ratio1'] >= 0.01) | (df1['ratio2'] >= 0.01))]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
list1.append(df1.shape[0])

df1 = pd.DataFrame(list5, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[(df1['RR'] > 1) & ((df1['ratio1'] >= 0.01) | (df1['ratio2'] >= 0.01))]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
list1.append(df1.shape[0])

df1 = pd.DataFrame(list6, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[(df1['RR'] > 1) & ((df1['ratio1'] >= 0.01) | (df1['ratio2'] >= 0.01))]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
list1.append(df1.shape[0])

df1 = pd.DataFrame(list7, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[(df1['RR'] > 1) & ((df1['ratio1'] >= 0.01) | (df1['ratio2'] >= 0.01))]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
list1.append(df1.shape[0])

df1 = pd.DataFrame(list7, columns=['disease1', 'disease2', 'description1', 'description2', 'C1', 'C2', 'C12', 'Csimultineously', 'ratio1', 'ratio2', 'RR', 'p_value', 'category1', 'category2'])
df1 = df1[df1['RR'] > 1]
df1 = df1[df1['p_value'] < (0.05/df1.shape[0])]
print(df1.shape[0])
list1.append(df1.shape[0])

list1.reverse()
plt.bar(range(0, len(list1)), list1, color='blue')
plt.xticks(range(0, len(list1)), ['original', '5 year', '4 year', '3 year', '2 year', '1 year', 'half year', '1 day'], fontsize=12)
# plt.show()
plt.savefig('a.pdf', bbox_inches='tight')