import pandas as pd
import xlrd


disease_class = dict()
disease_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
lis1 = df.values.tolist()
for each in lis1:
    code = each[0]
    c = each[2]
    disease_class[code] = c
    for each1 in code.split(';'):
        disease_merged[each1] = code


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

genetic_comorbidity = set()
df = pd.read_csv('../phenotype/comorbidity_filter.csv')
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    if (code1 in genetic_disease) & (code2 in genetic_disease):
        genetic_comorbidity.add((code1, code2))

within_category_comorbidity = set()
between_category_comorbidity = set()
for each in genetic_comorbidity:
    [code1, code2] = each
    c1 = disease_class[code1]
    c2 = disease_class[code2]
    if (c1 == 'Male genital organs') & (c2 == 'Male genital organs'):
        within_category_comorbidity.add((code1, code2))
    if (c1 != c2) & ((c1 == 'Male genital organs') | (c2 == 'Male genital organs')):
        between_category_comorbidity.add((code1, code2))

print('within category comorbidity: ' + str(len(within_category_comorbidity)))
print('cross category comorbidity: ' + str(len(between_category_comorbidity)))

gene_comorbidity = list()
with open('../overlap/comorbidity_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2, name1, name2 = str2[:4]
        c1 = disease_class[code1]
        c2 = disease_class[code2]
        set1 = set(str2[5].split(';'))
        if (code1, code2) in within_category_comorbidity:
            for each in set1:
                gene_comorbidity.append([each, code1, code2, name1, name2, c1, c2])
    infile.close()

pathway_comorbidity = list()
with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2, name1, name2 = str2[:4]
        c1 = disease_class[code1]
        c2 = disease_class[code2]
        set1 = set(str2[4].split(';'))
        if (code1, code2) in between_category_comorbidity:
            for each in set1:
                pathway_comorbidity.append([each, code1, code2, name1, name2, c1, c2])
    infile.close()


df = pd.DataFrame(gene_comorbidity, columns=['gene','code1', 'code2', 'name1', 'name2', 'c1', 'c2'])
df.to_csv('a.csv', index=False)
df = pd.DataFrame(pathway_comorbidity, columns=['pathway', 'code1', 'code2', 'name1', 'name2', 'c1', 'c2'])
df.to_csv('a1.csv', index=False)


for each in set(df['pathway']):
    df1 = df[df['pathway']==each]
    print([each, df1.shape[0]])