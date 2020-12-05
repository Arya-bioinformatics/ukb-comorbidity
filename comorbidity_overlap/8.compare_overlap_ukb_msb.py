from scipy import stats
import pandas as pd
import xlrd


if __name__ == '__main__':

    disease_merged = dict()
    df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
    list1 = df.values.tolist()
    for each in list1:
        for each1 in each[0].split(';'):
            disease_merged[each1] = each[0]

    genetic_disease = set()
    file = xlrd.open_workbook('../phenotype/geneAtlas_EMR.xlsx')
    sheet = file.sheets()[0]
    nrows = sheet.nrows
    for rownum in range(1, nrows):
        row = sheet.row_values(rownum)
        icd10 = str(row[0]).split('_')[-1]
        if icd10 in disease_merged:
            icd10 = disease_merged[icd10]
            genetic_disease.add(icd10)
    print(len(genetic_disease))

    genetic_comorbidity = set()
    df = pd.read_csv('../phenotype/comorbidity_filter.csv')
    list1 = df.values.tolist()
    for each in list1:
        code1, code2 = each[:2]
        if (code1 in genetic_disease) & (code2 in genetic_disease):
            genetic_comorbidity.add((code1, code2))


    icd9cm_icd10_map = dict()
    df = pd.read_table('../name_map/icd10_icd9cm_map.txt')
    list1 = df.values.tolist()
    for each in list1:
        icd10 = each[0]
        icd9cm = each[1]
        if len(icd10) > 3:
            icd10 = icd10[:3]
        set1 = set(icd9cm.split(';'))
        for each1 in set1:
            icd9cm_icd10_map[each1] = icd10

    msb_disease = set()
    msb_diseasePair = set()
    df = pd.read_excel('../name_map/msb_download/inline-supplementary-material-3.xls')
    list1 = df.values.tolist()
    for each in list1:
        icd9_code1 = each[0].replace('[', '').replace(']', '')
        icd9_code2 = each[2].replace('[', '').replace(']', '')
        if icd9_code1 not in icd9cm_icd10_map:
            continue
        if icd9_code2 not in icd9cm_icd10_map:
            continue
        icd10_code1 = icd9cm_icd10_map[icd9_code1]
        icd10_code2 = icd9cm_icd10_map[icd9_code2]
        if icd10_code1 in disease_merged:
            icd10_code1 = disease_merged[icd10_code1]
        if icd10_code2 in disease_merged:
            icd10_code2 = disease_merged[icd10_code2]
        msb_disease.add(icd10_code1)
        msb_disease.add(icd10_code2)
        if icd10_code1 < icd10_code2:
            msb_diseasePair.add((icd10_code1, icd10_code2))
        else:
            msb_diseasePair.add((icd10_code2, icd10_code1))


    ukb_comorbidity_snp = set()
    with open('../overlap/comorbidity_snp.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1, code2 = str2[:2]
            ukb_comorbidity_snp.add((code1, code2))
        infile.close()

    ukb_comorbidity_gene = set()
    with open('../overlap/comorbidity_gene.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1, code2 = str2[:2]
            ukb_comorbidity_gene.add((code1, code2))
        infile.close()

    ukb_comorbidity_ppi = set()
    with open('../overlap/comorbidity_ppi.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1, code2 = str2[:2]
            ukb_comorbidity_ppi.add((code1, code2))
        infile.close()

    ukb_comorbidity_pathway = set()
    with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1, code2 = str2[:2]
            ukb_comorbidity_pathway.add((code1, code2))
        infile.close()

    ukb_comorbidity_rg = set()
    with open('../overlap/comorbidity_rg.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1, code2 = str2[:2]
            ukb_comorbidity_rg.add((code1, code2))
        infile.close()

    print('------------- msb and ukb common used genetic diseases -----------------')
    set1 = genetic_disease & msb_disease
    print(len(set1))

    print('------------- ukb genetic and msb genetic overlap -----------------')
    set1 = genetic_comorbidity & msb_diseasePair
    print(len(set1))


    a = len(msb_diseasePair & (ukb_comorbidity_snp | ukb_comorbidity_gene | ukb_comorbidity_ppi | ukb_comorbidity_pathway | ukb_comorbidity_rg))
    b = len(msb_diseasePair & genetic_comorbidity) - a
    c = len((ukb_comorbidity_snp | ukb_comorbidity_gene | ukb_comorbidity_ppi | ukb_comorbidity_pathway | ukb_comorbidity_rg)) - a
    d = 8212 - a - b - c

    print(' -------- compare the overlap of msb and ukb ---------- ')
    [odds, p] = stats.fisher_exact([[a, b], [c, d]])
    print([odds, p])