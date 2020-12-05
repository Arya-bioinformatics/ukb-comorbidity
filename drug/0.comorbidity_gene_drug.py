import pandas as pd


drug_disease = dict()
with open('../drug_disease/ctd_drug_disease_info.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        drugid = str2[0]
        icd10 = str2[-2]
        evidence = str2[-1]
        # if evidence == '':
        #     continue
        if icd10 == '':
            continue
        if drugid not in drug_disease:
            drug_disease[drugid] = set()
        set1 = icd10.split(';')
        for each in set1:
            each1 = each.split('.')[0]
            drug_disease[drugid].add(each1)
    infile.close()


comorbidity_gene = dict()
df = pd.read_table('../overlap/comorbidity_gene.txt', dtype=str)
list1 = df.values.tolist()
for each in list1:
    code1 = each[0]
    code2 = each[1]
    set1 = set(each[4].split(';'))
    comorbidity_gene[(code1, code2)] = set1


gene_drugid = dict()
drug_id_name_map = dict()
df = pd.read_table('../drug_target/drugbank_target.txt', dtype=str)
list1 = df.values.tolist()
for each in list1:
    drugid = each[0]
    drugname = each[1]
    drug_id_name_map[drugid] = drugname
    gene_set = set(each[7].split(';'))
    for each1 in gene_set:
        if each1 not in gene_drugid:
            gene_drugid[each1] = set()
        gene_drugid[each1].add(drugid)


list1 = list()
for each in comorbidity_gene:
    code1, code2 = each
    set1 = comorbidity_gene[each]
    for each1 in set1:
        if each1 in gene_drugid:
            set2 = gene_drugid[each1]
            if len(set2) == 0:
                continue
            for each2 in set2:
                flag1 = 'N'
                flag2 = 'N'
                name = drug_id_name_map[each2]
                if each2 in drug_disease:
                    disease_set = drug_disease[each2]
                else:
                    disease_set = set()
                for each_disease in disease_set:
                    if len(each_disease) == 3:
                        if each_disease == code1:
                            flag1 = 'Y'
                        if each_disease == code2:
                            flag2 = 'Y'
                    elif len(each_disease) == 7:
                        [d1, d2] = each_disease.split('-')
                        if (code1 >= d1) & (code1 <= d2):
                            flag1 = 'Y'
                        if (code2 >= d1) & (code2 <= d2):
                            flag2 = 'Y'
                    if (flag1 == 'Y') & (flag2 == 'Y'):
                        break
                list1.append([code1, code2, each1, each2, name, ';'.join(disease_set), flag1, flag2])


print(len(list1))


df = pd.DataFrame(list1, columns=['code1', 'code2', 'gene', 'drug_id', 'drug_name', 'drug_indication', 'treat disease1', 'treat disease2'])

df.to_csv('a.csv', index=False)