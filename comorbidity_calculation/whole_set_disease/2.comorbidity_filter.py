import pandas as pd
from itertools import combinations

if __name__ == '__main__':

    disease_description = dict()
    disease_patient_number = dict()
    dict1 = dict()
    whole_set_disease = set()
    df = pd.read_excel('../phenotype/Disease_infomation.xlsx')
    list1 = df.values.tolist()
    for each in list1:
        code, name, num = each
        whole_set_disease.add(code)
        disease_description[code] = name
        disease_patient_number[code] = num
        list2 = list(code.split(';'))
        for each1 in list2:
            dict1[each1] = code


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
                if each[:3] in dict1:
                    each1 = dict1[each[:3]]
                else:
                    each1 = each
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

    print('all patient diagnosis: ' + str(len(patient_disease)))
    print('all patient date: ' + str(len(patient_diagnosis_date)))

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


    df = pd.read_csv('../phenotype/comorbidity_alldisease.csv')
    list1 = df.values.tolist()
    print('total disease pairs: ' + str(len(list1)))
    list2 = list()
    for each in list1:
        code1 = each[0]
        code2 = each[1]
        Ii = int(each[2])
        Ij = int(each[3])
        Cij = int(each[4])
        rr = float(each[5])
        p = float(each[6])
        name1 = disease_description[code1]
        name2 = disease_description[code2]
        key = code1 + '-' + code2
        if key in double_dis_pat:
            Cij_sim = len(double_dis_pat[key])
        else:
            Cij_sim = 0
        ratio1 = Cij_sim / Ii
        ratio2 = Cij_sim / Ij
        list2.append([code1, code2, name1, name2, Ii, Ij, Cij, Cij_sim, ratio1, ratio2, rr, p])

    df1 = pd.DataFrame(list2, columns=['disease1', 'disease2', 'description1', 'description2', 'Ii', 'Ij', 'Cij', 'Cij_sim', 'Cijsim/Ii', 'Cijsim/Ij', 'RR', 'p value'])


    df2 = df1.loc[(df1['RR'] > 1) & (df1['p value'] < (0.05 / df1.shape[0])) & ((df1['Cijsim/Ii'] > 0.01) | (df1['Cijsim/Ij'] > 0.01))]
    print(df2.shape[0])
    df2.to_csv('../phenotype/comorbidity_filter_alldisease.csv', index=False)


    df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
    subset_disease = set(df['ICD10'])
    common_disease = whole_set_disease & subset_disease
    print(' -------- commonly used disease ------------ ')
    print(len(common_disease))

    print(' ------- common disease comorbidity by whole set of diseases ------- ')
    df3 = df2[(df2['Ii']>410) & (df2['Ij']>410)]
    print(df3.shape[0])

    list1 = df3.values.tolist()
    set1 = set()
    for each in list1:
        code1, code2 = each[:2]
        if (code1 in common_disease) & (code2 in common_disease):
            set1.add((code1, code2))
    df4 = pd.read_csv('../phenotype/comorbidity_filter.csv')
    list1 = df4.values.tolist()
    set2 = set()
    for each in list1:
        code1, code2 = each[:2]
        if (code1 in common_disease) & (code2 in common_disease):
            set2.add((code1, code2))
    print('------- common disease comorbidity intersaction between using whole set of diseases and subset of diseases -------')
    print(len(set1))
    print(len(set2))
    print(len(set1 & set2))