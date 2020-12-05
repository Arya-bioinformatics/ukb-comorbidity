import pandas as pd
from scipy.stats import fisher_exact, ttest_ind
from matplotlib import pyplot as plt
import seaborn as sns
import os
import numpy as np


for flag in ['with', 'without']:
    print('\n\n' + flag + '\n\n')
    if flag == 'with':
        path1 = '../overlap/comorbidity_snp.txt'
        path2 = '../genome/disease_snp_sig_ld.txt'
        path3 = '../genome/varaint_info.txt'
        path4 = '../functional_annotation/all_variant.variant_function'
        path5 = 'a.csv'
        path6 = 'a.pdf'
        # continue
    if flag == 'without':
        path1 = '../overlap/comorbidity_snp_rmhla.txt'
        path2 = '../genome/disease_snp_sig_ld_rmhla.txt'
        path3 = '../genome/varaint_info_rmhla.txt'
        path4 = '../functional_annotation/all_variant_rmhla.variant_function'
        path5 = 'a1.csv'
        path6 = 'a1.pdf'
        # continue


    comorbidity_snp = set()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            comorbidity_snp |= set(str2[-1].split(';'))
        infile.close()

    disease_snp = set()
    with open(path2, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_snp |= set(str2[-1].split(';'))
        infile.close()

    all_snp = set()
    rsid_position_map = dict()
    with open(path3, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            all_snp.add(str2[0])
            if str2[3] < str2[4]:
                rsid_position_map[str2[1] + ':' + str2[2] + ':' + str2[3] + ':' + str2[4]] = str2[0]
            else:
                rsid_position_map[str2[1] + ':' + str2[2] + ':' + str2[4] + ':' + str2[3]] = str2[0]
        infile.close()

    # ----------------------- annovar --------------------- #
    annovar_gene_annotation = dict()
    with open(path4, 'r') as infile:
        for line in infile:
            str1 = line.strip('\r\n')
            str2 = str1.split()
            c = str2[0]
            rsid = str2[-1]
            if ('exonic' == c) | ('intronic' == c) | ('splicing' == c):
                annovar_gene_annotation[rsid] = 'Genic region'
            elif 'RNA' in c:
                annovar_gene_annotation[rsid] = 'Noncoding RNA'
            else:
                annovar_gene_annotation[rsid] = 'Intergenic region'
        infile.close()

    all_type = ['Noncoding RNA', 'Intergenic region', 'Genic region']
    nonassciated_type = dict()
    only_disease_type = dict()
    comorbidity_type = dict()
    count1 = 0
    count2 = 0
    count3 = 0
    for each in annovar_gene_annotation:
        type = annovar_gene_annotation[each]
        if each in comorbidity_snp:
            if type not in comorbidity_type:
                comorbidity_type[type] = 0
            comorbidity_type[type] += 1
            count1 += 1
        elif each in disease_snp:
            if type not in only_disease_type:
                only_disease_type[type] = 0
            only_disease_type[type] += 1
            count2 += 1
        else:
            if type not in nonassciated_type:
                nonassciated_type[type] = 0
            nonassciated_type[type] += 1
            count3 += 1

    list1 = list()
    for each in all_type:
        a = nonassciated_type[each]
        b = only_disease_type[each]
        c = comorbidity_type[each]
        list1.append([each, a, b, c, float(a) / count3, float(b) / count2, float(c) / count1])

    df = pd.DataFrame(list1, columns=['annotation', 'nonassociated', 'only_disease', 'Comorbidity',
                               'Ratio_of_nonassociated', 'Ratio_of_only_disease', 'Ratio_of_comorbidity'])
    df.to_csv(path5, index=False)

    # enrichment
    print('------------ compared to all -------------')
    for each in all_type:
        a = comorbidity_type[each]
        b = nonassciated_type[each] + only_disease_type[each]
        c = len(comorbidity_snp) - a
        d = len(all_snp) - a - b - c
        [odds, p] = fisher_exact([[a, b], [c, d]])
        print([each, odds, p])

    print('------------ compared to disease -------------')
    for each in all_type:
        a = comorbidity_type[each]
        b = only_disease_type[each]
        c = len(comorbidity_snp) - a
        d = len(disease_snp) - a - b - c
        [odds, p] = fisher_exact([[a, b], [c, d]])
        print([each, odds, p])



    # ----------------------- CADD --------------------- #
    variant_CADD_score = dict()
    with open('../functional_annotation/CADD_download/gnomad.genomes.r2.1.1.snv.tsv', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 2:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            if str2[2] < str2[3]:
                key = str2[0] + ':' + str2[1] + ':' + str2[2] + ':' + str2[3]
                if key in rsid_position_map:
                    rsid = rsid_position_map[key]
                    if rsid in variant_CADD_score:
                        print(rsid)
                    variant_CADD_score[rsid] = float(str2[5])
            else:
                key = str2[0] + ':' + str2[1] + ':' + str2[3] + ':' + str2[2]
                if key in rsid_position_map:
                    rsid = rsid_position_map[key]
                    if rsid in variant_CADD_score:
                        print(rsid)
                    variant_CADD_score[rsid] = float(str2[5])

        infile.close()

    print(len(variant_CADD_score))

    nonassociated_score = list()
    only_disease_score = list()
    comorbidity_score = list()
    for each in variant_CADD_score:
        score = variant_CADD_score[each]
        if each in comorbidity_snp:
            comorbidity_score.append(score)
        elif each in disease_snp:
            only_disease_score.append(score)
        else:
            nonassociated_score.append(score)

    # t test
    print('------------ compare to nonassociated --------------')
    [t, p] = ttest_ind(comorbidity_score, nonassociated_score)
    print([t, p])
    print('------------ compare to only disease --------------')
    [t, p] = ttest_ind(comorbidity_score, only_disease_score)
    print([t, p])
    print('------------ compare disease to nonassociated --------------')
    [t, p] = ttest_ind(only_disease_score, nonassociated_score)
    print([t, p])


    list1 = list()
    for each in nonassociated_score:
        list1.append(['Non-disease SNPs', each])
    for each in only_disease_score:
        list1.append(['Other disease-SNPs', each])
    for each in comorbidity_score:
        list1.append(['Comorbidity-SNPs', each])



    df = pd.DataFrame(list1, columns=['group', 'CADD score'])
    b = sns.violinplot(x="group", y="CADD score", data=df)
    b.set_xlabel("", fontsize=1)
    b.set_ylabel("CADD score", fontsize=16)
    b.tick_params(labelsize=16)
    b.set_xticklabels(['Non-disease SNPs', 'Other disease-SNPs', 'Comorbidity-SNPs'], rotation=45)
    plt.savefig(path6, bbox_inches='tight')
    plt.close()


    # ----------------------- dbscSNV_RF_SCORE  --------------------- #
    variant_dbscSNV_RF_SCORE = dict()
    file_list = os.listdir('../functional_annotation/dbscSNV_download/')
    for each in file_list:
        if 'chr' not in each:
            continue
        with open('../functional_annotation/dbscSNV_download/' + each, 'r') as infile:
            for i, line in enumerate(infile):
                if i < 1:
                    continue
                str1 = line.strip('\r\n')
                str2 = str1.split('\t')
                if str2[2] < str2[3]:
                    key = ':'.join(str2[:4])
                else:
                    key = ':'.join([str2[0], str2[1], str2[3], str2[2]])
                if key in rsid_position_map:
                    if str2[-1] == '.':
                        continue
                    variant_dbscSNV_RF_SCORE[rsid_position_map[key]] = float(str2[-1])
            infile.close()

    nonassociated_score = list()
    only_disease_score = list()
    comorbidity_score = list()
    for each in variant_dbscSNV_RF_SCORE:
        score = variant_dbscSNV_RF_SCORE[each]
        if each in comorbidity_snp:
            comorbidity_score.append(score)
        elif each in disease_snp:
            only_disease_score.append(score)
        else:
            nonassociated_score.append(score)

    # t test
    print('------------ compare to nonassociated --------------')
    [t, p] = ttest_ind(comorbidity_score, nonassociated_score)
    print([t, p])
    print([np.mean(comorbidity_score), np.mean(nonassociated_score)])
    print('------------ compare to only disease --------------')
    [t, p] = ttest_ind(comorbidity_score, only_disease_score)
    print([t, p])
    print([np.mean(comorbidity_score), np.mean(only_disease_score)])
    print('------------ compare disease to nonassociated --------------')
    [t, p] = ttest_ind(only_disease_score, nonassociated_score)
    print([t, p])
    print([np.mean(only_disease_score), np.mean(nonassociated_score)])


    list1 = list()
    for each in nonassociated_score:
        list1.append(['Non-disease SNPs', each])
    for each in only_disease_score:
        list1.append(['Other disease-SNPs', each])
    for each in comorbidity_score:
        list1.append(['Comorbid-SNPs', each])

    df = pd.DataFrame(list1, columns=['group', 'Splicing score'])
    b = sns.violinplot(x="group", y="Splicing score", data=df)
    b.set_xlabel("", fontsize=1)
    b.set_ylabel("Splicing score", fontsize=16)
    b.tick_params(labelsize=16)
    b.set_xticklabels(['Non-disease SNPs', 'Other disease-SNPs', 'Comorbid-SNPs'], rotation=45)
    plt.savefig(path6, bbox_inches='tight')
    plt.close()