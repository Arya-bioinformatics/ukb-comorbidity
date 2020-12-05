import pandas as pd

if __name__ == '__main__':

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


    all_ICD10 = set()
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
                all_ICD10.add(each1)
        infile.close()



    icd10_phecode = dict()
    df = pd.read_csv('../phenotype/phecode_icd10.csv', dtype=str)
    print(df.shape[0])
    df = df.dropna(axis=0, how='any')
    print(df.shape[0])
    list1 = df[['ICD10', 'PheCode']].values.tolist()
    for each in list1:
        icd10 = each[0]
        if len(icd10) != 3:
            continue
        if icd10[0] > 'N':
            continue
        phecode = str(each[1])
        icd10_phecode[icd10] = phecode

    list1 = list()
    all_phecode = set()
    for each in all_ICD10:
        if each not in icd10_phecode:
            print(each)
            list1.append([each, each])
            all_phecode.add(each)
        else:
            phecode = icd10_phecode[each]
            list1.append([each, phecode])
            all_phecode.add(phecode)
    df = pd.DataFrame(list1, columns=['icd10', 'phecode'])

    merged_icd10_name = list()
    dict1 = dict()
    for each in all_phecode:
        df1 = df[df['phecode'] == each]
        list1 = list(df1['icd10'])
        list1.sort()
        for each1 in list1:
            dict1[each1] = ';'.join(list1)
        list2 = list()
        for each1 in list1:
            name = disease_code_description[each1]
            list2.append(name)
        merged_icd10_name.append([';'.join(list1), ';'.join(list2)])

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


    with open('../phenotype/merged_icd10.txt', 'w+') as outfile:
        outfile.write('ICD10' + '\t' + 'Description' + '\t' + 'Number of patients' + '\n')
        for each in merged_icd10_name:
            patient_num = len(disease_patient[each[0]])
            outfile.write('\t'.join(each) + '\t' + str(patient_num) + '\n')
        outfile.close()