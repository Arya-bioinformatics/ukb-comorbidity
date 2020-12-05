import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats import multitest
import numpy as np
import xlrd

disease_merged = dict()
disease_category = dict()
all_category = set()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    for each1 in each[0].split(';'):
        disease_merged[each1] = each[0]
    disease_category[each[0]] = each[2]
    all_category.add(each[2])

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


list2 = list()
df = pd.read_csv('../phenotype/comorbidity_filter.csv')
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    if (code1 in genetic_disease) & (code2 in genetic_disease):
        list2.append([code1, code2, disease_category[code1], disease_category[code2]])
genetic_df = pd.DataFrame(list2, columns=['code1', 'code2', 'c1', 'c2'])

snp_interpreted = list()
df = pd.read_table('../overlap/comorbidity_snp.txt')
list1 = df.values.tolist()
for each in list1:
    snp_interpreted.append([each[0], each[1], disease_category[each[0]], disease_category[each[1]]])

gene_interpreted = list()
df = pd.read_table('../overlap/comorbidity_gene.txt')
list1 = df.values.tolist()
for each in list1:
    gene_interpreted.append([each[0], each[1], disease_category[each[0]], disease_category[each[1]]])

ppi_interpreted = list()
df = pd.read_table('../overlap/comorbidity_ppi.txt')
list1 = df.values.tolist()
for each in list1:
    ppi_interpreted.append([each[0], each[1], disease_category[each[0]], disease_category[each[1]]])

pathway_interpreted = list()
df = pd.read_table('../overlap/comorbidity_pathway.txt')
list1 = df.values.tolist()
for each in list1:
    pathway_interpreted.append([each[0], each[1], disease_category[each[0]], disease_category[each[1]]])

rg_interpreted = list()
df = pd.read_table('../overlap/comorbidity_rg.txt')
list1 = df.values.tolist()
for each in list1:
    rg_interpreted.append([each[0], each[1], disease_category[each[0]], disease_category[each[1]]])


snp_df = pd.DataFrame(snp_interpreted, columns=['code1', 'code2', 'c1', 'c2'])
gene_df = pd.DataFrame(gene_interpreted, columns=['code1', 'code2', 'c1', 'c2'])
ppi_df = pd.DataFrame(ppi_interpreted, columns=['code1', 'code2', 'c1', 'c2'])
pathway_df = pd.DataFrame(pathway_interpreted, columns=['code1', 'code2', 'c1', 'c2'])
rg_df = pd.DataFrame(rg_interpreted, columns=['code1', 'code2', 'c1', 'c2'])

inner = list()
outer = list()
p_list = list()
df1 = genetic_df[genetic_df['c1'] == genetic_df['c2']]
genetic_inner = df1.shape[0]
genetic_outer = genetic_df.shape[0] - genetic_inner

df1 = snp_df[snp_df['c1'] == snp_df['c2']]
a = df1.shape[0]
b = snp_df.shape[0] - a
inner.append(float(a)/genetic_inner)
outer.append(float(b)/genetic_outer)
c = genetic_inner - a
d = genetic_df.shape[0] - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]])
print([odds, p])
p_list.append(p)


df1 = gene_df.loc[gene_df['c1'] == gene_df['c2']]
a = df1.shape[0]
b = gene_df.shape[0] - a
inner.append(float(a) / genetic_inner)
outer.append(float(b) / genetic_outer)
c = genetic_inner - a
d = genetic_df.shape[0] - a - b - c
[odds, p] = fisher_exact([[a,b], [c,d]])
print([odds, p])
p_list.append(p)


df1 = ppi_df.loc[ppi_df['c1'] == ppi_df['c2']]
a = df1.shape[0]
b = ppi_df.shape[0] - a
inner.append(float(a) / genetic_inner)
outer.append(float(b) / genetic_outer)
c = genetic_inner - a
d = genetic_df.shape[0] - a - b - c
[odds, p] = fisher_exact([[a, b], [c, d]])
print([odds, p])
p_list.append(p)


df1 = pathway_df.loc[pathway_df['c1'] == pathway_df['c2']]
a = df1.shape[0]
b = pathway_df.shape[0] - a
inner.append(float(a) / genetic_inner)
outer.append(float(b) / genetic_outer)
c = genetic_inner - a
d = genetic_df.shape[0] - a - b - c
[odds, p] = fisher_exact([[a, b], [c, d]])
print([odds, p])
p_list.append(p)


df1 = rg_df.loc[rg_df['c1'] == rg_df['c2']]
a = df1.shape[0]
b = rg_df.shape[0] - a
inner.append(float(a) / genetic_inner)
outer.append(float(b) / genetic_outer)
c = genetic_inner - a
d = genetic_df.shape[0] - a - b - c
[odds, p] = fisher_exact([[a, b], [c, d]])
print([odds, p])
p_list.append(p)

q_list = multitest.fdrcorrection(p_list)[1]
print(q_list)
print(inner)
print(outer)



x = np.arange(5)
inner = [i*100 for i in inner]
outer = [i * 100 for i in outer]
bar_width = 0.30

fig = plt.figure(figsize=(6, 4))
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.5
plt.rcParams["font.family"] = "Arial"
plt.bar(x, inner, color="#F8766D", edgecolor='black', linewidth=0.5, width=0.3, label='within-category', hatch="///")
plt.bar(x+bar_width, outer, color="#F8766D", edgecolor='black', linewidth=0.5, width=0.3, label='cross-category', hatch="...")
plt.legend(fontsize=15, loc=2)
plt.xticks(x + bar_width / 2, ['SNP', 'Gene', 'PPI', 'Pathway', 'Genetic\ncorrelation'])
plt.ylabel("Ratio of comorbidities\nexplained by genetics, %", fontsize=15)
plt.tick_params(labelsize=15)
plt.ylim(0, 35)
plt.text(-0.2, 4, 'P=4.6e-4', fontsize=12)
plt.text(0.8, 19, 'P=0.024', fontsize=12)
plt.text(1.8, 24, 'P=4.6e-17', fontsize=12)
plt.text(2.8, 26, 'P=8.3e-17', fontsize=12)
plt.text(3.8, 25, 'P=0.677', fontsize=12)

plt.savefig('b.pdf', bbox_inches='tight')