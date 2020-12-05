import os

# with HLA
list1 = list()
file_list = os.listdir('../genome/magma_file/with_HLA/')
for each in file_list:
    if each.endswith('.txt'):
        key = each.replace('.txt', '').split('_')[-1]
        if len(key) > 3:
            print(key)
            continue
        if key[0] > 'N':
            print(key)
            continue
        str1 = '../magma/magma_v1.07b/magma ' \
               '--bfile ../magma/g1000_eur/g1000_eur ' \
               'synonyms=../magma/g1000_eur/g1000_eur.synonyms ' \
               '--gene-annot ../genome/magma_file/snp.genes.annot ' \
               '--pval ../genome/magma_file/with_HLA/' + each + \
               ' N=452264 duplicate=first --out ../genome/magma_file/annotate_with_HLA/' + key
        list1.append(str1)

list2 = list()
i = 0
for each in list1:
    list2.append(each)
    if len(list2) == 30:
        with open('../genome/magma_file/sh_gene_ananlysis/gene_analysis_with_HLA' + str(i) + '.sh', 'w+') as outfile:
            outfile.write('\n'.join(list2) + '\n')
            outfile.close()
        i += 1
        list2 = list()

with open('../genome/magma_file/sh_gene_ananlysis/gene_analysis_with_HLA' + str(i) + '.sh', 'w+') as outfile:
    outfile.write('\n'.join(list2) + '\n')
    outfile.close()

# without HLA
list1 = list()
file_list = os.listdir('../genome/magma_file/without_HLA/')
for each in file_list:
    if each.endswith('.txt'):
        key = each.replace('.txt', '').split('_')[-1]
        if len(key) > 3:
            print(key)
            continue
        if key[0] > 'N':
            print(key)
            continue
        str1 = '../magma/magma_v1.07b/magma ' \
               '--bfile ../magma/g1000_eur/g1000_eur ' \
               'synonyms=../magma/g1000_eur/g1000_eur.synonyms ' \
               '--gene-annot ../genome/magma_file/snp_rmhla.genes.annot ' \
               '--pval ../genome/magma_file/without_HLA/' + each + \
               ' N=452264 duplicate=first --out ../genome/magma_file/annotate_without_HLA/' + key
        list1.append(str1)

list2 = list()
i = 0
for each in list1:
    list2.append(each)
    if len(list2) == 30:
        with open('../genome/magma_file/sh_gene_ananlysis/gene_analysis_without_HLA' + str(i) + '.sh', 'w+') as outfile:
            outfile.write('\n'.join(list2) + '\n')
            outfile.close()
        i += 1
        list2 = list()

with open('../genome/magma_file/sh_gene_ananlysis/gene_analysis_without_HLA' + str(i) + '.sh', 'w+') as outfile:
    outfile.write('\n'.join(list2) + '\n')
    outfile.close()