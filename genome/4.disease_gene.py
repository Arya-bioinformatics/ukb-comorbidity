import pandas as pd
import os
from statsmodels.stats import multitest

for flag in ['with', 'without']:
    print(' ----------------- ' + flag + '----------------- ')
    if flag == 'with':
        path1 = '../genome/disease_snp_sig_ld.txt'
        path2 = '../genome/disease_gene_mapped.txt'
        path3 = '../genome/disease_gene_eqtl.txt'
        path4 = '../genome/magma_file/annotate_with_HLA/'
        path5 = '../genome/disease_gene_magma.txt'
        path6 = '../genome/disease_gene.txt'
    elif flag == 'without':
        path1 = '../genome/disease_snp_sig_ld_rmhla.txt'
        path2 = '../genome/disease_gene_mapped_rmhla.txt'
        path3 = '../genome/disease_gene_eqtl_rmhla.txt'
        path4 = '../genome/magma_file/annotate_without_HLA/'
        path5 = '../genome/disease_gene_magma_rmhla.txt'
        path6 = '../genome/disease_gene_rmhla.txt'

    # # ------------------ disease gene by genome coordination mapping ----------------- #
    disease_snp_ld = dict()
    all_disease_snp = set()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_snp_ld[str2[0]] = set(str2[1].split(';'))
            all_disease_snp |= set(str2[1].split(';'))
        infile.close()

    variant_position = dict()
    with open('../LD_Calculation/1000_genome_download/1000genomes_from_plink2_grch37/1000genomes.bim', 'r') as infile:
        for line in infile:
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            rsid = str2[1]
            if rsid not in all_disease_snp:
                continue
            chromesome = str2[0]
            position = str2[3]
            variant_position[rsid] = [chromesome, position]
        infile.close()

    with open('../genome/varaint_info.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            rsid = str2[0]
            if rsid not in all_disease_snp:
                continue
            chromesome = str2[1]
            position = str2[2]
            variant_position[rsid] = [chromesome, position]
        infile.close()

    # gene position
    gene_position = dict()
    all_gene = set()
    with open('../genome/grch37_download/genome_assemblies/gene_assembly_grch37.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\n')
            str2 = str1.split('\t')
            gene_id = str2[0]
            chr = str2[2]
            start = str2[3]
            end = str2[4]
            if (start == '*') | (end == '*'):
                continue
            type = str2[5]
            if type == 'pseudo':
                continue
            all_gene.add(gene_id)
            if chr in gene_position:
                gene_position[chr].append([start, end, gene_id])
            else:
                gene_position[chr] = list()
                gene_position[chr].append([start, end, gene_id])
        infile.close()

    snp_gene_mapping = dict()
    for each in all_disease_snp:
        if each in variant_position:
            [chro, pos] = variant_position[each]
            list1 = gene_position[chro]
            for s, e, gene in list1:
                if (int(pos) >= int(s) - 2000) & (int(pos) <= int(e) + 500):
                    if each not in snp_gene_mapping:
                        snp_gene_mapping[each] = set()
                    snp_gene_mapping[each].add(gene)

    disease_gene_mapped = dict()
    for each in disease_snp_ld:
        set1 = disease_snp_ld[each]
        set2 = set()
        for each1 in set1:
            if each1 in snp_gene_mapping:
                set2 |= snp_gene_mapping[each1]
        disease_gene_mapped[each] = set2
    print('total disease with gene by genome coordination mapping: ' + str(len(disease_gene_mapped)))
    with open(path2, 'w+') as outfile:
        for each in disease_gene_mapped:
            outfile.write(each + '\t' + ';'.join(disease_gene_mapped[each]) + '\n')
        outfile.close()

    # ------------------ disease gene by eqtl ----------------- #
    ensembl_gene_id = dict()
    with open('../eqtl/GTEX_download/custom_hgnc_ensembl.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            ncbiid = str2[4]
            ensemblid = str2[5]
            if ensemblid == '':
                continue
            if ncbiid == '':
                continue
            ensembl_gene_id[ensemblid] = ncbiid
        infile.close()

    snp_eqtl_gene_df = pd.read_table('../eqtl/GTEX_download/all_snp_gene_eqtl.txt')
    snp_eqtl_gene_df.columns = ['rsid', 'variant', 'gene', 'pvalue_fe']
    print('snp eqtl no correction: ' + str(snp_eqtl_gene_df.shape[0]))
    snp_eqtl_gene_df['q'] = multitest.fdrcorrection(snp_eqtl_gene_df['pvalue_fe'])[1]
    snp_eqtl_gene_df = snp_eqtl_gene_df[snp_eqtl_gene_df['q'] < 0.05]
    print('snp eqtl correction: ' + str(snp_eqtl_gene_df.shape[0]))

    disease_gene_eqtl = dict()
    temp = set()
    for each in disease_snp_ld:
        set1 = set(snp_eqtl_gene_df.loc[snp_eqtl_gene_df['rsid'].isin(disease_snp_ld[each]), 'gene'])
        set2 = set()
        for each1 in set1:
            ensembl = each1.split('.')[0]
            if ensembl in ensembl_gene_id:
                set2.add(ensembl_gene_id[ensembl])
        disease_gene_eqtl[each] = set2 & all_gene
        temp |= set2
    print('total disease gene by eqtl: ' + str(len(temp)))
    print('total disease with gene by eqtl: ' + str(len(disease_gene_eqtl)))
    with open(path3, 'w+') as outfile:
        for each in disease_gene_eqtl:
            outfile.write(each + '\t' + ';'.join(disease_gene_eqtl[each]) + '\n')
        outfile.close()

    # ------------------ disease gene by magma ----------------- #
    icd10_merged = dict()
    df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
    list1 = df.values.tolist()
    for each in list1:
        set1 = set(each[0].split(';'))
        for each1 in set1:
            icd10_merged[each1] = each[0]

    file_list = os.listdir(path4)

    list1 = list()
    for each in file_list:
        if each.endswith('.genes.out'):
            icd10 = each.replace('.genes.out', '')
            if icd10 not in icd10_merged:
                continue
            merged = icd10_merged[icd10]
            with open(path4 + each, 'r') as infile:
                for i, line in enumerate(infile):
                    if i < 1:
                        continue
                    str1 = line.strip('\r\n')
                    str2 = str1.split()
                    gene, p = str2[0], str2[-1]
                    list1.append([merged, gene, float(p)])
            infile.close()

    disease_gene_magma = dict()
    df = pd.DataFrame(list1, columns=['disease', 'gene', 'p'])
    df['q'] = multitest.fdrcorrection(df['p'])[1]
    df = df[df['q']<0.05]
    all_disease = set(df['disease'])
    for each in all_disease:
        disease_gene_magma[each] = set(df.loc[df['disease']==each, 'gene']) & all_gene

    with open(path5, 'w+') as outfile:
        for each in disease_gene_magma:
            outfile.write(each + '\t' + ';'.join(disease_gene_magma[each]) + '\n')
        outfile.close()

    disease_gene = dict()
    for each in disease_gene_mapped:
        disease_gene[each] = disease_gene_mapped[each]
    for each in disease_gene_eqtl:
        if each in disease_gene:
            disease_gene[each] |= disease_gene_eqtl[each]
        else:
            disease_gene[each] = disease_gene_eqtl[each]
    for each in disease_gene_magma:
        if each in disease_gene:
            disease_gene[each] |= disease_gene_magma[each]
        else:
            disease_gene[each] = disease_gene_magma[each]

    with open(path6, 'w+') as outfile:
        for each in disease_gene:
            outfile.write(each + '\t' + ';'.join(disease_gene[each]) + '\n')
        outfile.close()