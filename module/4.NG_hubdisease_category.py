import pandas as pd
from scipy.stats import fisher_exact
from itertools import combinations
from statsmodels.stats import multitest
import xlrd


disease_merged = dict()
class_disease = dict()
disease_name = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    disease_name[each[0]] = each[1]
    for each1 in each[0].split(';'):
        disease_merged[each1] = each[0]
    if each[2] not in class_disease:
        class_disease[each[2]] = set()
    class_disease[each[2]].add(each[0])

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

dict1 = dict()
for each in class_disease:
    dict1[each] = class_disease[each] & genetic_disease
class_disease = dict1

df = pd.read_csv('../module/NG_network_module_Q0.3200458.csv')
NG_module_disease = dict()
list1 = df.values.tolist()
for each in list1:
    disease, module = each[0], each[3]
    if module not in NG_module_disease:
        NG_module_disease[module] = set()
    NG_module_disease[module].add(disease)

NG_module_Encategory = dict()
list1 = list()
for each in NG_module_disease:
    set1 = NG_module_disease[each]
    if len(set1) < 5:
        continue
    for each1 in class_disease:
        set2 = class_disease[each1]
        a = len(set1 & set2)
        b = len(set2) - a
        c = len(set1) - a
        d = len(genetic_disease) - a - b - c
        [odds, p] = fisher_exact([[a,b], [c,d]])
        list1.append([each, each1, odds, p])
df = pd.DataFrame(list1, columns=['module', 'class', 'odds', 'p'])
df = df[df['odds']>1]
df['q'] = multitest.fdrcorrection(df['p'])[1]
df = df[df['q']<0.05]
list1 = df.values.tolist()
for each in list1:
    module, category = each[:2]
    if module not in NG_module_Encategory:
        NG_module_Encategory[module] = set()
    NG_module_Encategory[module].add(category)

NG_module_hubDisease = dict()
NG_module_hubPPI = dict()
NG_module_hubPathway = dict()
comorbidity_set = set()
comorbidity_ppi = dict()
comorbidity_pathway = dict()
df = pd.read_table('../overlap/comorbidity_ppi.txt')
list1 = df.values.tolist()
for each in list1:
    comorbidity_set.add((each[0], each[1]))
    comorbidity_ppi[(each[0], each[1])] = set(each[5].split(';'))
df = pd.read_table('../overlap/comorbidity_pathway.txt')
list1 = df.values.tolist()
for each in list1:
    comorbidity_set.add((each[0], each[1]))
    comorbidity_pathway[(each[0], each[1])] = set(each[4].split(';'))

for each in NG_module_disease:
    set1 = NG_module_disease[each]
    if len(set1) < 5:
        continue
    module_diseasePair = set()
    for code1, code2 in combinations(set1, 2):
        if code1 < code2:
            module_diseasePair.add((code1, code2))
        else:
            module_diseasePair.add((code2, code1))
    module_comorbidity = module_diseasePair & comorbidity_set
    disease_comorbidity = dict()
    for each1 in set1:
        set2 = set()
        for each2 in module_comorbidity:
            if each1 in each2:
                set2.add(each2)
        disease_comorbidity[each1] = set2
    list1 = list()
    for each1 in disease_comorbidity:
        list1.append([each1, len(disease_comorbidity[each1])])
    df = pd.DataFrame(list1, columns=['disease', 'comorbidity number'])
    df = df.sort_values(by='comorbidity number', ascending=False)
    list1 = list(set(df['comorbidity number']))
    list1.sort(reverse=True)
    top_3_num = list1[:3]
    hub_disease = list()
    set2 = set()
    for each1 in df.values.tolist():
        disease, number = each1
        if number in top_3_num[:1]:
            hub_disease.append(disease)
            set2 |= disease_comorbidity[disease]
    if len(set2) / len(module_comorbidity) < 0.5:
        hub_disease = list()
        set2 = set()
        for each1 in df.values.tolist():
            disease, number = each1
            if number in top_3_num[:2]:
                hub_disease.append(disease)
                set2 |= disease_comorbidity[disease]
    if len(set2) / len(module_comorbidity) < 0.5:
        hub_disease = list()
        set2 = set()
        for each1 in df.values.tolist():
            disease, number = each1
            if number in top_3_num:
                hub_disease.append(disease)
                set2 |= disease_comorbidity[disease]
    if len(set2) / len(module_comorbidity) < 0.5:
        continue
    NG_module_hubDisease[each] = hub_disease
    hub_disease_ppi = list()
    hub_disease_pathway = list()
    for each1 in hub_disease:
        ppi_list = list()
        for each2 in disease_comorbidity[each1]:
            if each2 not in comorbidity_ppi:
                continue
            ppi_list += list(comorbidity_ppi[each2])
        list1 = list()
        for each2 in set(ppi_list):
            list1.append([each2, ppi_list.count(each2)])
        df = pd.DataFrame(list1, columns=['ppi', 'count'])
        df = df.sort_values(by='count', ascending=False)
        list1 = list(set(df['count']))
        list1.sort(reverse=True)
        top_5_num = list1[:5]
        top_ppi = list()
        for each2 in df.values.tolist():
            if each2[1] in top_5_num:
                top_ppi.append(each2[0])
        hub_disease_ppi.append(disease_name[each1] + ': ' + ', '.join(top_ppi))

        pathway_list = list()
        for each2 in disease_comorbidity[each1]:
            if each2 not in comorbidity_pathway:
                continue
            pathway_list += list(comorbidity_pathway[each2])
        list1 = list()
        for each2 in set(pathway_list):
            list1.append([each2, pathway_list.count(each2)])
        df = pd.DataFrame(list1, columns=['pathway', 'count'])
        list1 = list(set(df['count']))
        list1.sort(reverse=True)
        top_5_num = list1[:5]
        top_pathway = list()
        for each2 in df.values.tolist():
            if each2[1] in top_5_num:
                top_pathway.append(each2[0])
        hub_disease_pathway.append(disease_name[each1] + ': ' + ', '.join(top_pathway))
    NG_module_hubPPI[each] = '\n'.join(hub_disease_ppi)
    NG_module_hubPathway[each] = '\n'.join(hub_disease_pathway)

list1 = list()
module_order = list(NG_module_disease.keys())
module_order.sort()
for each in module_order:
    module = 'NG-module' + str(each)
    size = len(NG_module_disease[each])
    disease = ' | '.join(NG_module_disease[each])
    if each in NG_module_Encategory:
        temp = list(NG_module_Encategory[each])
        temp.sort()
        category = '-'.join(temp)
    else:
        category = '/'
    if each in NG_module_hubDisease:
        hubs = ' | '.join(NG_module_hubDisease[each])
    else:
        hubs = '/'
    if each in NG_module_hubPPI:
        ppi = NG_module_hubPPI[each]
    else:
        ppi = '/'
    if each in NG_module_hubPathway:
        pathway = NG_module_hubPathway[each]
    else:
        pathway = '/'
    list1.append([module, size, disease, category, hubs, ppi, pathway])
df = pd.DataFrame(list1, columns=['Comorbidity Module', 'Module size', 'Disease', 'Featured categories', 'Hub diseases', 'Top 5 PPI shared by hub diseases', 'Top 5 Pathway shared by hub diseases'])
df.to_csv('a1.csv', index=False)