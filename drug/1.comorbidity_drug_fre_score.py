import pandas as pd
import xlrd

if __name__ == '__main__':

    drug_code_name = dict()
    df = pd.read_table('/../phenotype/coding4.tsv', dtype=str)
    list1 = df.values.tolist()
    for each in list1:
        code = each[0]
        name = each[1]
        drug_code_name[code] = name

    disease_category = dict()
    file = xlrd.open_workbook('/../phenotype/disease_class_manual.xlsx')
    sheet = file.sheets()[0]
    nrows = sheet.nrows
    for rownum in range(0, nrows):
        row = sheet.row_values(rownum)
        icd10 = str(row[0])
        category = str(row[2])
        disease_category[icd10] = category

    double_disease_patient = dict()
    df = pd.read_table('/../phenotype/double_disease_patient_phecode.txt', dtype=str)
    list1 = df.values.tolist()
    for each in list1:
        double_disease_patient[each[0]] = set(each[1].split(';'))

    patient_drug = dict()
    with open('/../phenotype/field20003.csv', 'r') as infile:
        lineNum = 1
        for line in infile:
            if lineNum == 1:
                lineNum = 2
                continue
            str1 = line.strip('\r\n')
            str1 = str1.replace('"', '')
            str2 = str1.split(',')
            set1 = set(str2[1:])
            if '' in set1:
                set1.remove('')
            if '99999' in set1:
                set1.remove('99999')
            if len(set1) == 0:
                continue
            patient_drug[str2[0]] = set1
        infile.close()

    df = pd.read_csv('/../phenotype/comorbidity_filter.csv', dtype=str)
    df1 = df[['disease1', 'disease2']]
    comorbidity_list = df1.values.tolist()

    list1 = list()
    for each in comorbidity_list:
        code1 = each[0]
        code2 = each[1]
        c1 = disease_category[code1]
        c2 = disease_category[code2]
        patient_set = double_disease_patient[code1 + '-' + code2]
        for each1 in patient_set:
            if each1 in patient_drug:
                drug_set = patient_drug[each1]
                for each2 in drug_set:
                    list1.append([code1, code2, c1, c2, each1, each2])

    df = pd.DataFrame(list1, columns=['disease1', 'disease2', 'c1', 'c2', 'patient', 'drug'])

    comorbidity_drug_score = list()
    i = 0
    for each in comorbidity_list:
        print(i)
        i = i + 1
        code1 = each[0]
        code2 = each[1]
        c1 = disease_category[code1]
        c2 = disease_category[code2]
        drug_set = set(df.loc[(df['disease1'] == code1) & (df['disease2'] == code2), 'drug'])

        temp1 = df.loc[(df['disease1'] == code1) & (df['disease2'] == code2), ]
        b = set(temp1['patient'])
        if len(b) == 0:
            continue

        for each1 in drug_set:
            if each1 not in drug_code_name:
                continue
            name = drug_code_name[each1]
            temp = temp1.loc[temp1['drug'] == each1, 'patient']
            a = set(temp)
            x1 = len(a)
            x2 = len(b)
            comorbidity_drug_score.append([code1, code2, name, str(float(x1) / x2)])
    df2 = pd.DataFrame(comorbidity_drug_score, columns=['disease1', 'disease2', 'drug', 'frequency'])
    df2.to_csv('/../overlap/drug/comorbidity_drug_score.csv', index=False)