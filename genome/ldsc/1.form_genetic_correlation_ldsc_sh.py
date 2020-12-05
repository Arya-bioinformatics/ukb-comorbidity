import os
import xlrd


# 读入prevalence
disease_prevalence = dict()
file = xlrd.open_workbook('../geneAtlas/41588_2018_248_MOESM3_ESM.xlsx')
sheet = file.sheets()[0]
nrows = sheet.nrows
ncols = sheet.ncols
for rownum in range(4, nrows):
    row = sheet.row_values(rownum)
    key = str(row[0])
    prevalence = str(row[6])
    disease_prevalence[key] = prevalence


# with HLA
has_done = set()
file_list = os.listdir('../genome/ldsc_file/ldsc_rg_with_HLA/')
for each in file_list:
    if os.path.getsize('../genome/ldsc_file/ldsc_rg_with_HLA/' + each):
        has_done.add(each.replace('.log', ''))

root_dir = '../genome/ldsc_file/sumstats_with_HLA/'
all_file = os.listdir(root_dir)
sumstats_file = list()
for each in all_file:
    if each.endswith('sumstats.gz'):
        sumstats_file.append(each)

list1 = list()
for i in range(0, len(sumstats_file) - 1):
    for j in range(i, len(sumstats_file)):
        key1 = sumstats_file[i].replace('.sumstats.gz', '')
        key2 = sumstats_file[j].replace('.sumstats.gz', '')
        if (key1 + '~' + key2) in has_done:
            continue
        prev1 = disease_prevalence[key1]
        prev2 = disease_prevalence[key2]
        str1 = 'python ../LDSC/ldsc-master/ldsc.py' \
               ' --rg ../genome/ldsc_file/sumstats_with_HLA/' + sumstats_file[i] + ',' + '../genome/ldsc_file/sumstats_with_HLA/' + sumstats_file[j] + \
               ' --ref-ld-chr ../LDSC/eur_w_ld_chr/' \
               ' --w-ld-chr ../LDSC/eur_w_ld_chr/' \
               ' --out ../genome/ldsc_file/ldsc_rg_with_HLA/' + key1 + '~' + key2 + ' --samp-prev ' + prev1 + ',' + prev2 + \
               ' --pop-prev ' + prev1 + ',' + prev2
        list1.append(str1)

list2 = list()
i = 0
for each in list1:
    list2.append(each)
    if len(list2) == 5000:
        with open('../genome/ldsc_file/sh_rg/ldsc_with_HLA' + str(i) +'.sh', 'w+') as outfile:
            outfile.write('\n'.join(list2) + '\n')
            outfile.close()
        list2 = list()
        i += 1

with open('../genome/ldsc_file/sh_rg/ldsc_with_HLA' + str(i) +'.sh', 'w+') as outfile:
    outfile.write('\n'.join(list2) + '\n')
    outfile.close()



# without HLA
has_done = set()
file_list = os.listdir('../genome/ldsc_file/ldsc_rg_without_HLA/')
for each in file_list:
    if os.path.getsize('../genome/ldsc_file/ldsc_rg_without_HLA/' + each):
        has_done.add(each.replace('.log', ''))

root_dir = '../genome/ldsc_file/sumstats_without_HLA/'
all_file = os.listdir(root_dir)
sumstats_file = list()
for each in all_file:
    if each.endswith('sumstats.gz'):
        sumstats_file.append(each)

list1 = list()
for i in range(0, len(sumstats_file) - 1):
    for j in range(i, len(sumstats_file)):
        key1 = sumstats_file[i].replace('.sumstats.gz', '')
        key2 = sumstats_file[j].replace('.sumstats.gz', '')
        if (key1 + '~' + key2) in has_done:
            continue
        prev1 = disease_prevalence[key1]
        prev2 = disease_prevalence[key2]
        str1 = 'python ../LDSC/ldsc-master/ldsc.py' \
               ' --rg ../genome/ldsc_file/sumstats_without_HLA/' + sumstats_file[i] + ',' + '../genome/ldsc_file/sumstats_without_HLA/' + sumstats_file[j] + \
               ' --ref-ld-chr ../LDSC/eur_w_ld_chr/' \
               ' --w-ld-chr ../LDSC/eur_w_ld_chr/' \
               ' --out ../genome/ldsc_file/ldsc_rg_without_HLA/' + key1 + '~' + key2 + ' --samp-prev ' + prev1 + ',' + prev2 + \
               ' --pop-prev ' + prev1 + ',' + prev2
        list1.append(str1)

list2 = list()
i = 0
for each in list1:
    list2.append(each)
    if len(list2) == 5000:
        with open('../genome/ldsc_file/sh_rg/ldsc_without_HLA' + str(i) +'.sh', 'w+') as outfile:
            outfile.write('\n'.join(list2) + '\n')
            outfile.close()
        list2 = list()
        i += 1

with open('../genome/ldsc_file/sh_rg/ldsc_without_HLA' + str(i) +'.sh', 'w+') as outfile:
    outfile.write('\n'.join(list2) + '\n')
    outfile.close()