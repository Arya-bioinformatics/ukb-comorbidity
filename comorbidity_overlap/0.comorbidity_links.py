import pandas as pd
from statsmodels.stats import multitest
import os


disease_merged = dict()
disease_name = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    for each1 in each[0].split(';'):
        disease_merged[each1] = each[0]
    disease_name[each[0]] = each[1]

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


for flag in ['with', 'without']:
    if flag == 'with':
        path1 = '../genome/disease_snp_sig_ld.txt'
        path2 = '../genome/disease_gene.txt'
        path3 = '../genome/disease_ppi.txt'
        path4 = '../genome/disease_pathway.txt'
        path5 = '../genome/ldsc_file/diseasePair_rg.txt'
        path6 = '../overlap/comorbidity_snp.txt'
        path7 = '../overlap/comorbidity_gene.txt'
        path8 = '../overlap/comorbidity_ppi.txt'
        path9 = '../overlap/comorbidity_pathway.txt'
        path10 = '../overlap/comorbidity_rg.txt'
        path11 = '../overlap/comorbidity_overlap.txt'
        path12 = '../genome/filter_file/with_HLA/'
    elif flag == 'without':
        path1 = '../genome/disease_snp_sig_ld_rmhla.txt'
        path2 = '../genome/disease_gene_rmhla.txt'
        path3 = '../genome/disease_ppi_rmhla.txt'
        path4 = '../genome/disease_pathway_rmhla.txt'
        path5 = '../genome/ldsc_file/diseasePair_rg_rmhla.txt'
        path6 = '../overlap/comorbidity_snp_rmhla.txt'
        path7 = '../overlap/comorbidity_gene_rmhla.txt'
        path8 = '../overlap/comorbidity_ppi_rmhla.txt'
        path9 = '../overlap/comorbidity_pathway_rmhla.txt'
        path10 = '../overlap/comorbidity_rg_rmhla.txt'
        path11 = '../overlap/comorbidity_overlap_rmhla.txt'
        path12 = '../genome/filter_file/without_HLA/'

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


    df = pd.read_csv('../phenotype/comorbidity_filter.csv')
    comorbidity_snp = dict()
    comorbidity_gene = dict()
    comorbidity_ppi = dict()
    comorbidity_pathway = dict()
    comorbidity_rg = dict()
    list1 = df.values.tolist()
    list2 = list()
    for each in list1:
        code1, code2 = each[:2]
        des1 = disease_name[code1]
        des2 = disease_name[code2]
        if (code1 in disease_snp) & (code2 in disease_snp):
            set1 = disease_snp[code1] & disease_snp[code2]
            if len(set1) != 0:
                comorbidity_snp[(code1, code2, des1, des2)] = set1
        if (code1 in disease_gene) & (code2 in disease_gene):
            set1 = disease_gene[code1] & disease_gene[code2]
            if len(set1) != 0:
                comorbidity_gene[(code1, code2, des1, des2)] = set1
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
                comorbidity_ppi[(code1, code2, des1, des2)] = set4
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
                if (code1, code2, des1, des2) not in comorbidity_pathway:
                    comorbidity_pathway[(code1, code2, des1, des2)] = set()
                comorbidity_pathway[(code1, code2, des1, des2)] |= set3
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
                    if (code1, code2, des1, des2) not in comorbidity_pathway:
                        comorbidity_pathway[(code1, code2, des1, des2)] = set()
                    comorbidity_pathway[(code1, code2, des1, des2)] |= set3

        if (code1 + '~' + code2) in diseasePair_rg:
            [rg, p] = diseasePair_rg[code1 + '~' + code2]
            comorbidity_rg[(code1, code2, des1, des2)] = [rg, p]

    list1 = list()
    for each in comorbidity_rg:
        list1.append([each[0], each[1], each[2], each[3], comorbidity_rg[each][0], comorbidity_rg[each][1]])
    df1 = pd.DataFrame(list1, columns=['code1', 'code2', 'des1', 'des2', 'rg', 'p'])
    df1['q'] = multitest.fdrcorrection(df1['p'])[1]
    df1 = df1[df1['q'] < 0.05]
    list1 = df1.values.tolist()
    comorbidity_rg = dict()
    for each in list1:
        comorbidity_rg[(each[0], each[1], each[2], each[3])] = [str(each[4]), str(each[5]), str(each[6])]

    with open(path6, 'w+') as outfile:
        outfile.write('code1\tcode2\tdescription1\tdescription2\tsnp\n')
        for each in comorbidity_snp:
            outfile.write('\t'.join(each) + '\t' + ';'.join(comorbidity_snp[each]) + '\n')
        outfile.close()

    with open(path7, 'w+') as outfile:
        outfile.write('code1\tcode2\tdescription1\tdescription2\tgene id\tgene symbol\n')
        for each in comorbidity_gene:
            list1 = list(comorbidity_gene[each])
            list1.sort()
            list2 = list()
            for each1 in list1:
                list2.append(gene_id_symbol_map[each1])
            outfile.write('\t'.join(each) + '\t' + ';'.join(list1) + '\t' + ';'.join(list2) + '\n')
        outfile.close()

    with open(path8, 'w+') as outfile:
        outfile.write('code1\tcode2\tdescription1\tdescription2\tppi(gene id)\tppi(gene symbol)\n')
        for each in comorbidity_ppi:
            list1 = list(comorbidity_ppi[each])
            list1.sort()
            list2 = list()
            for each1 in list1:
                gene1, gene2 = each1.split('~')
                list2.append(gene_id_symbol_map[gene1] + '~' + gene_id_symbol_map[gene2])
            outfile.write('\t'.join(each) + '\t' + ';'.join(list1) + '\t' + ';'.join(list2) + '\n')
        outfile.close()

    with open(path9, 'w+') as outfile:
        outfile.write('code1\tcode2\tdescription1\tdescription2\tpathway\n')
        for each in comorbidity_pathway:
            outfile.write('\t'.join(each) + '\t' + ';'.join(comorbidity_pathway[each]) + '\n')
        outfile.close()

    with open(path10, 'w+') as outfile:
        outfile.write('code1\tcode2\tdescription1\tdescription2\trg\tp\tq\n')
        for each in comorbidity_rg:
            outfile.write('\t'.join(each) + '\t' + '\t'.join(comorbidity_rg[each]) + '\n')
        outfile.close()

    file_list = os.listdir(path12)
    genetic_disease = set()
    for each in file_list:
        icd10 = each.split('_')[-2]
        if icd10 in disease_merged:
            genetic_disease.add(disease_merged[icd10])
    print('genetic disease count: ' + str(len(genetic_disease)))

    genetic_comorbidity = set()
    list1 = df.values.tolist()
    for each in list1:
        code1, code2 = each[:2]
        if (code1 in genetic_disease) & (code2 in genetic_disease):
            genetic_comorbidity.add(code1 + '-' + code2)
    print('genetic comorbidity count: ' + str(len(genetic_comorbidity)))