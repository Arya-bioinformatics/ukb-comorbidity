import pandas as pd
from statsmodels.stats import multitest
import os
from itertools import combinations



disease_merged = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    for each1 in each[0].split(';'):
        disease_merged[each1] = each[0]

gene_id_symbol_map = dict()
with open('../genome/grch37_download/genome_assemblies/gene_assembly_grch37.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        gene_id_symbol_map[str2[0]] = str2[1]
    infile.close()

pathway_gene = dict()
with open('../genome/MSigDB/c2.cp.v6.2.entrez.gmt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        pathway = str.lower(str2[0])
        gene = set(str2[2:])
        pathway_gene[pathway] = gene
    infile.close()

comorbidity_set = set()
df = pd.read_csv('../phenotype/comorbidity_filter.csv')
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    comorbidity_set.add(code1 + '-' + code2)


for flag in ['with', 'without']:
    if flag == 'with':
        path1 = '../genome/disease_snp_sig_ld.txt'
        path2 = '../genome/disease_gene.txt'
        path3 = '../genome/disease_ppi.txt'
        path4 = '../genome/disease_pathway.txt'
        path5 = '../genome/ldsc_file/diseasePair_rg.txt'
        path6 = '../overlap/noncomorbidity_snp.txt'
        path7 = '../overlap/noncomorbidity_gene.txt'
        path8 = '../overlap/noncomorbidity_ppi.txt'
        path9 = '../overlap/noncomorbidity_pathway.txt'
        path10 = '../overlap/noncomorbidity_rg.txt'
        path11 = '../genome/filter_file/with_HLA/'
    elif flag == 'without':
        path1 = '../genome/disease_snp_sig_ld_rmhla.txt'
        path2 = '../genome/disease_gene_rmhla.txt'
        path3 = '../genome/disease_ppi_rmhla.txt'
        path4 = '../genome/disease_pathway_rmhla.txt'
        path5 = '../genome/ldsc_file/diseasePair_rg_rmhla.txt'
        path6 = '../overlap/noncomorbidity_snp_rmhla.txt'
        path7 = '../overlap/noncomorbidity_gene_rmhla.txt'
        path8 = '../overlap/noncomorbidity_ppi_rmhla.txt'
        path9 = '../overlap/noncomorbidity_pathway_rmhla.txt'
        path10 = '../overlap/noncomorbidity_rg_rmhla.txt'
        path11 = '../genome/filter_file/without_HLA/'

    file_list = os.listdir(path11)
    genetic_disease = set()
    for each in file_list:
        icd10 = each.split('_')[-2]
        if icd10 in disease_merged:
            genetic_disease.add(disease_merged[icd10])
    print('genetic disease count: ' + str(len(genetic_disease)))

    disease_snp = dict()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_snp[str2[0]] = set(str2[1].split(';'))
        infile.close()

    disease_gene = dict()
    with open(path2, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_gene[str2[0]] = set(str2[1].split(';'))
        infile.close()

    disease_ppi = dict()
    with open(path3, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_ppi[str2[0]] = set(str2[1].split(';'))
        infile.close()

    disease_pathway = dict()
    with open(path4, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_pathway[str2[0]] = set(str2[1].split(';'))
        infile.close()

    diseasePair_rg = dict()
    with open(path5, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1, code2 = str2[0].split('~')
            rg, p = float(str2[1]), float(str2[2])
            if code1 in disease_merged:
                code1 = disease_merged[code1]
            if code2 in disease_merged:
                code2 = disease_merged[code2]
            if code1 < code2:
                key = code1 + '~' + code2
            else:
                key = code2 + '~' + code1
            if key in diseasePair_rg:
                [rg1, p1] = diseasePair_rg[key]
                if p1 > p:
                    diseasePair_rg[key] = [rg, p]
            else:
                diseasePair_rg[key] = [rg, p]
        infile.close()

    noncomorbidity_snp = dict()
    noncomorbidity_gene = dict()
    noncomorbidity_ppi = dict()
    noncomorbidity_pathway = dict()
    noncomorbidity_rg = dict()
    for code1, code2 in combinations(genetic_disease, 2):
        if ((code1, code2) in comorbidity_set) | ((code2, code1) in comorbidity_set):
            continue
        if (code1 in disease_snp) & (code2 in disease_snp):
            set1 = disease_snp[code1] & disease_snp[code2]
            if len(set1) != 0:
                if code1 < code2:
                    noncomorbidity_snp[(code1, code2)] = set1
                else:
                    noncomorbidity_snp[(code2, code1)] = set1
        if (code1 in disease_gene) & (code2 in disease_gene):
            set1 = disease_gene[code1] & disease_gene[code2]
            if len(set1) != 0:
                if code1 < code2:
                    noncomorbidity_gene[(code1, code2)] = set1
                else:
                    noncomorbidity_gene[(code2, code1)] = set1
        if (code1 in disease_ppi) & (code2 in disease_ppi):
            set1 = disease_gene[code1]
            set2 = disease_gene[code2]
            set3 = disease_ppi[code1] & disease_ppi[code2]
            set4 = set()
            for each1 in set3:
                gene1, gene2 = each1.split('~')
                if (len(set([gene1, gene2]) & set1) == 2) | (len(set([gene1, gene2]) & set2) == 2):
                    continue
                if ((gene1 in set1) & (gene2 in set2)) | ((gene2 in set1) & (gene1 in set2)):
                    set4.add(each1)
            if len(set4) != 0:
                if code1 < code2:
                    noncomorbidity_ppi[(code1, code2)] = set4
                else:
                    noncomorbidity_ppi[(code2, code1)] = set4
        if (code1 in disease_pathway) & (code2 in disease_gene):
            set1 = disease_pathway[code1]
            set2 = disease_gene[code2]
            set3 = set()
            for each1 in set1:
                geneSet = pathway_gene[each1]
                intersection = geneSet & set2
                if len(intersection) != 0:
                    set3.add(each1)
            if len(set3) != 0:
                if code1 < code2:
                    if (code1, code2) not in noncomorbidity_pathway:
                        noncomorbidity_pathway[(code1, code2)] = set()
                    noncomorbidity_pathway[(code1, code2)] |= set3
                else:
                    if (code2, code1) not in noncomorbidity_pathway:
                        noncomorbidity_pathway[(code2, code1)] = set()
                    noncomorbidity_pathway[(code2, code1)] |= set3
        if (code1 in disease_gene) & (code2 in disease_pathway):
            set1 = disease_pathway[code2]
            set2 = disease_gene[code1]
            set3 = set()
            for each1 in set1:
                geneSet = pathway_gene[each1]
                intersection = geneSet & set2
                if len(intersection) != 0:
                    set3.add(each1)
            if len(set3) != 0:
                if code1 < code2:
                    if (code1, code2) not in noncomorbidity_pathway:
                        noncomorbidity_pathway[(code1, code2)] = set()
                    noncomorbidity_pathway[(code1, code2)] |= set3
                else:
                    if (code2, code1) not in noncomorbidity_pathway:
                        noncomorbidity_pathway[(code2, code1)] = set()
                    noncomorbidity_pathway[(code2, code1)] |= set3
        if (code1 + '~' + code2) in diseasePair_rg:
            [rg, p] = diseasePair_rg[code1 + '~' + code2]
            noncomorbidity_rg[(code1, code2)] = [rg, p]
        elif (code2 + '~' + code1) in diseasePair_rg:
            [rg, p] = diseasePair_rg[code2 + '~' + code1]
            noncomorbidity_rg[(code2, code1)] = [rg, p]

    list1 = list()
    for each in noncomorbidity_rg:
        list1.append([each[0], each[1], noncomorbidity_rg[each][0], noncomorbidity_rg[each][1]])
    df1 = pd.DataFrame(list1, columns=['code1', 'code2', 'rg', 'p'])
    df1['q'] = multitest.fdrcorrection(df1['p'])[1]
    df1 = df1[df1['q'] < 0.05]
    list1 = df1.values.tolist()
    noncomorbidity_rg = dict()
    for each in list1:
        noncomorbidity_rg[(each[0], each[1])] = [str(each[2]), str(each[3]), str(each[4])]

    with open(path6, 'w+') as outfile:
        outfile.write('code1\tcode2\tsnp\n')
        for each in noncomorbidity_snp:
            outfile.write('\t'.join(each) + '\t' + ';'.join(noncomorbidity_snp[each]) + '\n')
        outfile.close()

    with open(path7, 'w+') as outfile:
        outfile.write('code1\tcode2\tgene id\tgene symbol\n')
        for each in noncomorbidity_gene:
            list1 = list(noncomorbidity_gene[each])
            list1.sort()
            list2 = list()
            for each1 in list1:
                list2.append(gene_id_symbol_map[each1])
            outfile.write('\t'.join(each) + '\t' + ';'.join(list1) + '\t' + ';'.join(list2) + '\n')
        outfile.close()

    with open(path8, 'w+') as outfile:
        outfile.write('code1\tcode2\tppi(gene id)\tppi(gene symbol)\n')
        for each in noncomorbidity_ppi:
            list1 = list(noncomorbidity_ppi[each])
            list1.sort()
            list2 = list()
            for each1 in list1:
                gene1, gene2 = each1.split('~')
                list2.append(gene_id_symbol_map[gene1] + '~' + gene_id_symbol_map[gene2])
            outfile.write('\t'.join(each) + '\t' + ';'.join(list1) + '\t' + ';'.join(list2) + '\n')
        outfile.close()

    with open(path9, 'w+') as outfile:
        outfile.write('code1\tcode2\tpathway\n')
        for each in noncomorbidity_pathway:
            outfile.write('\t'.join(each) + '\t' + ';'.join(noncomorbidity_pathway[each]) + '\n')
        outfile.close()

    with open(path10, 'w+') as outfile:
        outfile.write('code1\tcode2\trg\tp\tq\n')
        for each in noncomorbidity_rg:
            outfile.write('\t'.join(each) + '\t' + '\t'.join(noncomorbidity_rg[each]) + '\n')
        outfile.close()