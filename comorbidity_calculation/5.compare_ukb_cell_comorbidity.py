import pandas as pd
from itertools import product
from matplotlib_venn import venn2, venn2_circles
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


if __name__ == '__main__':

    double_disease_patient = dict()
    with open('../phenotype/double_disease_patient_phecode.txt', 'r') as infile:
        lineNum = 1
        for line in infile:
            if lineNum == 1:
                lineNum = 2
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            double_disease_patient[str2[0]] = len(set(str2[1].split(';')))
        infile.close()

    ukb_diseasePair = dict()
    df = pd.read_csv('../phenotype/comorbidity.csv')
    print(df.shape[0])
    list1 = df.values.tolist()
    for each in list1:
        code1 = each[0]
        code2 = each[1]
        rr = each[-2]
        p = each[-1]
        c1 = each[2]
        c2 = each[3]
        c12 = each[4]
        key = code1 + '-' + code2
        if key not in double_disease_patient:
            ctime = 0
        else:
            ctime = double_disease_patient[key]
        x1 = float(ctime) / c1
        x2 = float(ctime) / c2
        ukb_diseasePair[(code1, code2)] = [rr, p, c1, c2, c12, ctime, x1, x2]

    df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
    ukb_disease = set(df['ICD10'])
    print('ukb icd10 disease: ' + str(len(ukb_disease)))
    list1 = df.values.tolist()
    disease_name = dict()
    for each in list1:
        code = each[0]
        name = each[1]
        disease_name[code] = name

    disease_merge = dict()
    df = pd.read_excel('../phenotype/ICD10_merged.xlsx')
    list1 = df.values.tolist()
    for each in list1:
        disease_merge[each[0]] = each[1]


    complex_disease = dict()
    complex = set()
    df = pd.read_excel('../name_map/cell_download/1-s2.0-S0092867413010246-mmc2.xls',encoding='ISO-8859-1')
    list1 = df.values.tolist()
    for each in list1:
        name = str.lower(each[0])
        str1 = each[2].replace(' ', '')
        set1 = set(str1.split(','))
        set2 = set()
        for each1 in set1:
            if len(each1) > 3:
                each2 = each1[:3]
            else:
                each2 = each1
            if each2 in disease_merge:
                each3 = disease_merge[each2]
            else:
                each3 = each2
            set2.add(each3)
        complex_disease[name] = set2
        complex |= set2


    mendelian_disease = dict()
    mendelian = set()
    df = pd.read_excel('../name_map/cell_download/1-s2.0-S0092867413010246-mmc3.xls',encoding='ISO-8859-1')
    list1 = df.values.tolist()
    for each in list1:
        name = str.lower(each[0])
        str1 = each[2].replace(' ', '')
        set1 = set(str1.split(','))
        set2 = set()
        for each1 in set1:
            if len(each1) > 3:
                each2 = each1[:3]
            else:
                each2 = each1
            if each2 in disease_merge:
                each3 = disease_merge[each2]
            else:
                each3 = each2
            set2.add(each3)
        mendelian_disease[name] = set2
        mendelian |= set2


    print('cell complex disease: ' + str(len(complex_disease)))
    print('cell mendelian disease: ' + str(len(mendelian_disease)))

    disease_overlap = ukb_disease & (complex | mendelian)
    print(complex & mendelian)
    disease_overlap = disease_overlap - (complex & mendelian)
    print('disease overlap: ' + str(len(disease_overlap)))
    print('complex disease overlap: ' + str(len(disease_overlap & complex)))
    print('mendelian disease overlap: ' + str(len(disease_overlap & mendelian)))


    cell_comorbidity = set()
    df = pd.read_excel('../name_map/cell_download/1-s2.0-S0092867413010246-mmc4.xls')
    list1 = df.values.tolist()
    for each in list1[1:]:
        disease1 = str.lower(each[0])
        disease2 = str.lower(each[1])
        if (disease1 in complex_disease) & (disease2 in mendelian_disease):
            set1 = complex_disease[disease1]
            set2 = mendelian_disease[disease2]
            set1 = set1 & disease_overlap
            set2 = set2 & disease_overlap
            list1 = list(product(set1, set2))
            for each1 in list1:
                [code1, code2] = each1
                if code1 < code2:
                    [rr, p, c1, c2, c12, ctime, x1, x2] = ukb_diseasePair[(code1, code2)]
                    if (x1<0.01) & (x2<0.01):
                        continue
                    cell_comorbidity.add((code1, code2))
                else:
                    [rr, p, c1, c2, c12, ctime, x1, x2] = ukb_diseasePair[(code2, code1)]
                    if (x1<0.01) & (x2<0.01):
                        continue
                    cell_comorbidity.add((code2, code1))
        elif (disease2 in complex_disease) & (disease1 in mendelian_disease):
            set1 = complex_disease[disease2]
            set2 = mendelian_disease[disease1]
            set1 = set1 & disease_overlap
            set2 = set2 & disease_overlap
            list1 = list(product(set1, set2))
            for each1 in list1:
                [code1, code2] = each1
                if code1 < code2:
                    [rr, p, c1, c2, c12, ctime, x1, x2] = ukb_diseasePair[(code1, code2)]
                    if (x1 < 0.01) & (x2 < 0.01):
                        continue
                    cell_comorbidity.add((code1, code2))
                else:
                    [rr, p, c1, c2, c12, ctime, x1, x2] = ukb_diseasePair[(code2, code1)]
                    if (x1 < 0.01) & (x2 < 0.01):
                        continue
                    cell_comorbidity.add((code2, code1))

    print('involved cell comorbidity: ' + str(len(cell_comorbidity)))

    involved_ukb_comorbidity = set()
    df = pd.read_csv('../phenotype/comorbidity_filter.csv')
    print(df.shape[0])
    list1 = df.values.tolist()
    for each in list1:
        disease1 = each[0]
        disease2 = each[1]
        rr = each[-4]
        p = each[-3]
        if disease1 not in disease_overlap:
            continue
        if disease2 not in disease_overlap:
            continue
        if (disease1 in complex) & (disease2 in mendelian):
            involved_ukb_comorbidity.add((disease1, disease2))
        if (disease2 in complex) & (disease1 in mendelian):
            involved_ukb_comorbidity.add((disease1, disease2))


    print('involved ukb comorbidity: ' + str(len(involved_ukb_comorbidity)))
    set1 = involved_ukb_comorbidity & cell_comorbidity
    print('comorbidity overlap: ' + str(len(set1)))

    temp = cell_comorbidity - involved_ukb_comorbidity
    list1 = list()
    for each in temp:
        [code1, code2] = each
        [rr, p, c1, c2, c12, ctime, x1, x2] = ukb_diseasePair[(code1, code2)]
        list1.append([code1, code2, rr, p, c1, c2, c12, ctime, x1, x2])

    df = pd.DataFrame(list1, columns=['code1', 'code2', 'rr', 'p', 'c1', 'c2', 'c12', 'ctime', 'ctime/c1', 'ctime/c2'])
    df.to_csv('a.csv', index=False)

    set1 = involved_ukb_comorbidity - cell_comorbidity
    for each in set1:
        code1 = each[0]
        code2 = each[1]
        name1 = disease_name[code1]
        name2 = disease_name[code2]
        print(code1 + '\t' + code2 + '\t' + name1 + '\t' + name2)
    print(len(set1))

    fig = plt.figure(figsize=(10, 6), dpi=300)
    c = venn2([involved_ukb_comorbidity, cell_comorbidity], ('UKB', 'Cell'), alpha=0.8)
    c.get_patch_by_id('10').set_color('#F8766D')
    c.get_patch_by_id('01').set_color('#00BA38')
    c.get_patch_by_id('10').set_edgecolor('black')
    c.get_patch_by_id('01').set_edgecolor('black')
    for text in c.set_labels:
        text.set_fontsize(30)
    for text in c.subset_labels:
        text.set_fontsize(30)
    plt.savefig("b.pdf", bbox_inches="tight")