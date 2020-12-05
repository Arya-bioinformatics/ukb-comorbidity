from scipy.stats import fisher_exact
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from itertools import combinations
import xlrd


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
        genetic_comorbidity.add(code1 + '-' + code2)

genetic_diseasePair = set()
for code1, code2 in combinations(genetic_disease, 2):
    if code1 < code2:
        genetic_diseasePair.add(code1 + '-' + code2)
    else:
        genetic_diseasePair.add(code2 + '-' + code1)

comorbidity_snp = set()
with open('../overlap/comorbidity_snp.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_snp.add(code1 + '-' + code2)
    infile.close()

comorbidity_gene = set()
with open('../overlap/comorbidity_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_gene.add(code1 + '-' + code2)
    infile.close()

comorbidity_ppi = set()
with open('../overlap/comorbidity_ppi.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_ppi.add(code1 + '-' + code2)
    infile.close()

comorbidity_pathway = set()
with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_pathway.add(code1 + '-' + code2)
    infile.close()

comorbidity_rg = set()
with open('../overlap/comorbidity_rg.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_rg.add(code1 + '-' + code2)
    infile.close()


noncomorbidity_snp = set()
with open('../overlap/noncomorbidity_snp.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        noncomorbidity_snp.add(code1 + '-' + code2)
    infile.close()

noncomorbidity_gene = set()
with open('../overlap/noncomorbidity_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        noncomorbidity_gene.add(code1 + '-' + code2)
    infile.close()

noncomorbidity_ppi = set()
with open('../overlap/noncomorbidity_ppi.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        noncomorbidity_ppi.add(code1 + '-' + code2)
    infile.close()

noncomorbidity_pathway = set()
with open('../overlap/noncomorbidity_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        noncomorbidity_pathway.add(code1 + '-' + code2)
    infile.close()

noncomorbidity_rg = set()
with open('../overlap/noncomorbidity_rg.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        noncomorbidity_rg.add(code1 + '-' + code2)
    infile.close()


ratio_list_1 = list()
ratio_list_2 = list()
significance_list = list()

# snp test
a = len(comorbidity_snp)
b = len(noncomorbidity_snp)
c = len(genetic_comorbidity) - a
d = len(genetic_diseasePair) - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
ratio_list_1.append(a/(a+c))
ratio_list_2.append(b/(b+d))
significance_list.append(p)

# gene test
a = len(comorbidity_gene)
b = len(noncomorbidity_gene)
c = len(genetic_comorbidity) - a
d = len(genetic_diseasePair) - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
ratio_list_1.append(a/(a+c))
ratio_list_2.append(b/(b+d))
significance_list.append(p)

# ppi test
a = len(comorbidity_ppi)
b = len(noncomorbidity_ppi)
c = len(genetic_comorbidity) - a
d = len(genetic_diseasePair) - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
ratio_list_1.append(a/(a+c))
ratio_list_2.append(b/(b+d))
significance_list.append(p)

# pathway test
a = len(comorbidity_pathway)
b = len(noncomorbidity_pathway)
c = len(genetic_comorbidity) - a
d = len(genetic_diseasePair) - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
ratio_list_1.append(a/(a+c))
ratio_list_2.append(b/(b+d))
significance_list.append(p)

# genetic correlation test
a = len(comorbidity_rg)
b = len(noncomorbidity_rg)
c = len(genetic_comorbidity) - a
d = len(genetic_diseasePair) - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
ratio_list_1.append(a/(a+c))
ratio_list_2.append(b/(b+d))
significance_list.append(p)

# total test
a = len(comorbidity_snp | comorbidity_gene | comorbidity_ppi | comorbidity_pathway | comorbidity_rg)
print(a)
b = len(noncomorbidity_snp | noncomorbidity_gene | noncomorbidity_ppi | noncomorbidity_pathway | noncomorbidity_rg)
c = len(genetic_comorbidity) - a
d = len(genetic_diseasePair) - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
ratio_list_1.append(a/(a+c))
ratio_list_2.append(b/(b+d))
significance_list.append(p)

print(ratio_list_1)
print(ratio_list_2)
print(significance_list)


mpl.rcParams["font.sans-serif"] = ["SimHei"]
mpl.rcParams["axes.unicode_minus"] = False

x = np.arange(6)
y = list()
y1 = list()
for each in ratio_list_1:
    y.append(each*100)
for each in ratio_list_2:
    y1.append(each*100)

bar_width = 0.30
tick_label = ["SNP", "Gene", "PPI", "Pathway", "Genetic\ncorrelation", 'Total']

fig = plt.figure(figsize=(6,4))
plt.rcParams["font.family"] = "Arial"
plt.bar(x, y, bar_width, align="center", color="#F8766D", edgecolor='black', linewidth=0.5, label="comorbidity")
plt.bar(x + bar_width, y1, bar_width, color="#00BA38", edgecolor='black', linewidth=0.5, align="center", label="non-comorbidity")
plt.ylim((0, 65))

plt.ylabel("Ratio of comorbidities\nexplained by genetics, %", fontsize=15)
plt.xticks(x + bar_width / 2, tick_label, fontsize=15, rotation=90)

plt.text(-0.42, 6, 'P=2.1e-20', fontsize=12)
plt.text(0.6, 20, 'P=2.6e-128', fontsize=12)
plt.text(1.6, 27, 'P=2.5e-98', fontsize=12)
plt.text(2.8, 27, 'P=8.2e-62', fontsize=12)
plt.text(3.9, 27, 'P=0', fontsize=12)
plt.text(4.5, 48, 'P=8.3e-271', fontsize=12)

plt.tick_params(labelsize=15)

plt.legend(fontsize=15, loc='upper left')
plt.savefig('a.pdf', bbox_inches='tight')