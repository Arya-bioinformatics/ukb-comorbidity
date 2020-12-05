import os
import pandas as pd


icd10_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    set1 = set(each[0].split(';'))
    for each1 in set1:
        icd10_merged[each1] = each[0]


df = pd.read_csv('../phenotype/comorbidity_filter.csv')
all_comor_disease = set(df['disease1']) | set(df['disease2'])

# with HLA
root_dir = '../genome/filter_file/with_HLA/'
file_list = os.listdir(root_dir)

genetic_disease = set()
dict1 = dict()
for each in file_list:
    icd10 = each.split('_')[-2]
    if icd10 in icd10_merged:
        genetic_disease.add(icd10_merged[icd10])
        dict1[icd10] = each

print('genetic disease count: ' + str(len(genetic_disease)))

all_variant_count = 0
with open('../genome/varaint_info.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        all_variant_count += 1
    infile.close()

print('all variant count: ' + str(all_variant_count))


disease_snp_sig = dict()
j = 0
for each in all_comor_disease:
    if each not in genetic_disease:
        continue
    print(j)
    j = j + 1
    list1 = each.split(';')
    for each1 in list1:
        if each1 not in dict1:
            continue
        file = dict1[each1]
        with open(root_dir + file, 'r') as infile:
            for i, line in enumerate(infile):
                if i < 1:
                    continue
                str1 = line.strip('\r\n')
                str2 = str1.split('\t')
                rsid = str2[0]
                p_value = str2[4]
                if float(p_value) < (0.05/(all_variant_count * len(genetic_disease))):
                    if each not in disease_snp_sig:
                        disease_snp_sig[each] = set()
                    disease_snp_sig[each].add(rsid)
            infile.close()


with open('../genome/disease_snp_sig.txt', 'w+') as outfile:
    outfile.write('icd10_code' + '\t' + 'snp' + '\n')
    for each in disease_snp_sig:
        outfile.write(each + '\t' + ';'.join(disease_snp_sig[each]) + '\n')
    outfile.close()


# without HLA
root_dir = '../genome/filter_file/without_HLA/'
file_list = os.listdir(root_dir)

genetic_disease = set()
dict1 = dict()
for each in file_list:
    icd10 = each.split('_')[-2]
    if icd10 in icd10_merged:
        genetic_disease.add(icd10_merged[icd10])
        dict1[icd10] = each

print('genetic disease count: ' + str(len(genetic_disease)))

all_variant_count = 0
with open('../genome/varaint_info_rmhla.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        all_variant_count += 1
    infile.close()

print('all variant count: ' + str(all_variant_count))


disease_snp_sig = dict()
j = 0
for each in all_comor_disease:
    if each not in genetic_disease:
        continue
    print(j)
    j = j + 1
    list1 = each.split(';')
    for each1 in list1:
        if each1 not in dict1:
            continue
        file = dict1[each1]
        with open(root_dir + file, 'r') as infile:
            for i, line in enumerate(infile):
                if i < 1:
                    continue
                str1 = line.strip('\r\n')
                str2 = str1.split('\t')
                rsid = str2[0]
                p_value = str2[4]
                if float(p_value) < (0.05/(all_variant_count * len(genetic_disease))):
                    if each not in disease_snp_sig:
                        disease_snp_sig[each] = set()
                    disease_snp_sig[each].add(rsid)
            infile.close()


with open('../genome/disease_snp_sig_rmhla.txt', 'w+') as outfile:
    outfile.write('icd10_code' + '\t' + 'snp' + '\n')
    for each in disease_snp_sig:
        outfile.write(each + '\t' + ';'.join(disease_snp_sig[each]) + '\n')
    outfile.close()