import os
import xlrd


# phenotype case
phenotype_case = dict()
file = xlrd.open_workbook('../geneAtlas/41588_2018_248_MOESM3_ESM.xlsx')
sheet = file.sheets()[0]
nrows = sheet.nrows
for rownum in range(4, nrows):
    row = sheet.row_values(rownum)
    key = str(row[0])
    case = str(row[5]).split('.')[0]
    phenotype_case[key] = case

root_dir = '../genome/ldsc_file/with_HLA/'
file_list = os.listdir(root_dir)
list1 = list()
for each in file_list:
    key = each.replace('.txt', '')
    case = phenotype_case[key]
    total = '452264'
    control = str(int(total) - int(case))
    str1 = 'python ../LDSC/ldsc-master/munge_sumstats.py' \
           ' --sumstats ../genome/ldsc_file/with_HLA/' + each + \
           ' --N ' + total + \
           ' --N-cas ' + case + \
           ' --N-con ' + control + \
           ' --out ../genome/ldsc_file/sumstats_with_HLA/' + each.replace('.txt', '') + \
           ' --chunksize 500000' \
           ' --merge-alleles ../LDSC/w_hm3.snplist'
    list1.append(str1)

with open('../genome/ldsc_file/sh_munge_sumstats/munge_sumstats_with_HLA.sh', 'w+') as outfile:
    outfile.write('\n'.join(list1))
    outfile.close()

os.system('cd ../genome/ldsc_file/sh_munge_sumstats/\nsplit -l 100 munge_sumstats_with_HLA.sh -d munge_sumstats_with_HLA\n')


# remove hLA
root_dir = '../genome/ldsc_file/without_HLA/'
file_list = os.listdir(root_dir)
list1 = list()
for each in file_list:
    key = each.replace('.txt', '')
    case = phenotype_case[key]
    total = '452264'
    control = str(int(total) - int(case))
    str1 = 'python ../LDSC/ldsc-master/munge_sumstats.py' \
           ' --sumstats ../genome/ldsc_file/without_HLA/' + each + \
           ' --N ' + total + \
           ' --N-cas ' + case + \
           ' --N-con ' + control + \
           ' --out ../genome/ldsc_file/sumstats_without_HLA/' + each.replace('.txt', '') + \
           ' --chunksize 500000' \
           ' --merge-alleles ../LDSC/w_hm3.snplist'
    list1.append(str1)


with open('../genome/ldsc_file/sh_munge_sumstats/munge_sumstats_without_HLA.sh', 'w+') as outfile:
    outfile.write('\n'.join(list1))
    outfile.close()

os.system('cd ../genome/ldsc_file/sh_munge_sumstats/\nsplit -l 100 munge_sumstats_without_HLA.sh -d munge_sumstats_without_HLA\n')