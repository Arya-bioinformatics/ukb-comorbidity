import os
import gzip
import xlrd


all_variant = set()
variant_allele = dict()
with open('../genome/varaint_info.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        snp = str2[0]
        all_variant.add(snp)
        variant_allele[snp] = [str2[3], str2[4]]
    infile.close()

# remove hla
all_variant_rmhla = set()
variant_allele_rmhla = dict()
with open('../genome/varaint_info_rmhla.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        snp = str2[0]
        all_variant_rmhla.add(snp)
        variant_allele_rmhla[snp] = [str2[3], str2[4]]
    infile.close()

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


root_dir = '../geneAtlas/gwas_download/my_used/results/'
file_list = os.listdir(root_dir)
file_list.sort()

i = 0
for each in file_list[:50]:
    if each.startswith('.'):
        continue
    else:
        print(i)
        i = i + 1
        print(each)
        case = phenotype_case[each]
        each_phenotype_file = root_dir + each
        all_chr_summary_file = os.listdir(each_phenotype_file)
        list1 = list()
        list2 = list()
        list3 = list()
        list4 = list()
        list5 = list()
        list6 = list()
        list1.append('snp' + '\t' + 'allele' + '\t' + 'beta' + '\t' + 'se' + '\t' + 'p_value')
        list2.append('snp' + '\t' + 'A1' + '\t' + 'A2' + '\t' + 'sample_size' + '\t' + 'p_value' + '\t' + 'beta')
        list3.append('SNP' + '\t' + 'P')
        list4.append('snp' + '\t' + 'allele' + '\t' + 'beta' + '\t' + 'se' + '\t' + 'p_value')
        list5.append('snp' + '\t' + 'A1' + '\t' + 'A2' + '\t' + 'sample_size' + '\t' + 'p_value' + '\t' + 'beta')
        list6.append('SNP' + '\t' + 'P')
        for each_chr_file in all_chr_summary_file:
            if each_chr_file.endswith('.csv.gz'):
                chr_file = gzip.open(each_phenotype_file + '/' + each_chr_file)
                content = str(chr_file.read(), encoding="utf8")
                all_line = content.split('\n')
                for line in all_line[1:]:
                    str1 = line.strip('\r\n')
                    str2 = str1.split(' ')
                    rsid = str2[0]
                    if rsid in all_variant:
                        if len(str2) == 5:
                            A1 = str2[1]
                            beta = str2[2]
                            p = str2[4]
                            A2 = (set(variant_allele[rsid]) - set(A1)).pop()
                            list1.append('\t'.join(str2))
                            list2.append('\t'.join([rsid, A1, A2, case, p, beta]))
                            list3.append('\t'.join([rsid, p]))
                        else:
                            A1 = str2[1]
                            beta = str2[3]
                            se = str2[4]
                            p = str2[5]
                            A2 = (set(variant_allele[rsid]) - set(A1)).pop()
                            list1.append('\t'.join([rsid, A1, beta, se, p]))
                            list2.append('\t'.join([rsid, A1, A2, case, p, beta]))
                            list3.append('\t'.join([rsid, p]))
                    if rsid in all_variant_rmhla:
                        if len(str2) == 5:
                            A1 = str2[1]
                            beta = str2[2]
                            p = str2[4]
                            A2 = (set(variant_allele_rmhla[rsid]) - set(A1)).pop()
                            list4.append('\t'.join(str2))
                            list5.append('\t'.join([rsid, A1, A2, case, p, beta]))
                            list6.append('\t'.join([rsid, p]))
                        else:
                            A1 = str2[1]
                            beta = str2[3]
                            se = str2[4]
                            p = str2[5]
                            A2 = (set(variant_allele_rmhla[rsid]) - set(A1)).pop()
                            list4.append('\t'.join([rsid, A1, beta, se, p]))
                            list5.append('\t'.join([rsid, A1, A2, case, p, beta]))
                            list6.append('\t'.join([rsid, p]))
        with open('../genome/filter_file/with_HLA/' + each + '_filter.txt', 'w+') as outfile:
            outfile.write('\n'.join(list1))
            outfile.close()
        with open('../genome/ldsc_file/with_HLA/' + each + '.txt', 'w+') as outfile:
            outfile.write('\n'.join(list2))
            outfile.close()
        with open('../genome/magma_file/with_HLA/' + each + '.txt', 'w+') as outfile:
            outfile.write('\n'.join(list3))
            outfile.close()
        with open('../genome/filter_file/without_HLA/' + each + '_filter.txt', 'w+') as outfile:
            outfile.write('\n'.join(list4))
            outfile.close()
        with open('../genome/ldsc_file/without_HLA/' + each + '.txt', 'w+') as outfile:
            outfile.write('\n'.join(list5))
            outfile.close()
        with open('../genome/magma_file/without_HLA/' + each + '.txt', 'w+') as outfile:
            outfile.write('\n'.join(list6))
            outfile.close()