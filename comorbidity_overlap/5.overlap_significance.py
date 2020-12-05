from scipy.stats import fisher_exact
import pandas as pd
from statsmodels.stats import multitest
import random

def pathwayEnrich(disease_gene_dict):
    disease_pathway = set()
    with open('../kegg_pathway/download/code_make/pathway_category.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            category = str2[0]
            if category != 'Human Diseases':
                continue
            name = str.lower(str2[2]).replace(' ', '_')
            if name == 'leishmaniasis':
                name = 'leishmania_infection'
            if name == 'graft-versus-host_disease':
                name = 'graft_versus_host_disease'
            disease_pathway.add('kegg_' + name)
        infile.close()

    pathway_geneid = dict()
    all_gene = set()
    with open('../genome/MSigDB/c2.cp.v6.2.entrez.gmt', 'r') as infile:
        for line in infile:
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            pathway = str.lower(str2[0])
            if pathway in disease_pathway:
                continue
            geneid = set(str2[2:])
            if len(geneid) > 200:
                continue
            pathway_geneid[pathway] = geneid
            all_gene = all_gene | geneid
        infile.close()

    all_pathway = list(pathway_geneid.keys())
    list1 = list()
    for each in disease_gene_dict:
        disease_gene_set = disease_gene_dict[each]
        for each1 in all_pathway:
            if len(disease_gene_set) < 3:
                continue
            pathway_gene_set = pathway_geneid[each1]
            a = len(disease_gene_set & pathway_gene_set)
            b = len(pathway_gene_set) - a
            c = len(disease_gene_set & all_gene) - a
            d = len(all_gene) - a - b - c
            [odds, p] = fisher_exact([[a, b], [c, d]], alternative='greater')
            list1.append([each, each1, odds, p])
    return list1, pathway_geneid


all_variant_count = 0
with open('../genome/varaint_info.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        all_variant_count += 1
    infile.close()

all_gene_count = 0
all_gene = set()
with open('../genome/grch37_download/genome_assemblies/gene_assembly_grch37.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        if str2[-1] == 'pseudo':
            continue
        all_gene_count += 1
        all_gene.add(str2[0])
    infile.close()

ppi_dict = dict()
with open('../genome/PPI_biogrid/PPI_physical.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        gene1, gene2 = str1.split('\t')
        if gene1 not in ppi_dict:
            ppi_dict[gene1] = set()
        ppi_dict[gene1].add(gene2)
        if gene2 not in ppi_dict:
            ppi_dict[gene2] = set()
        ppi_dict[gene2].add(gene1)
    infile.close()


disease_snp = dict()
with open('../genome/disease_snp_sig_ld.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease_snp[str2[0]] = set(str2[1].split(';'))
    infile.close()

disease_gene = dict()
all_disease_genes = set()
with open('../genome/disease_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease_gene[str2[0]] = set(str2[1].split(';'))
        all_disease_genes |= set(str2[1].split(';'))
    infile.close()

disease_ppi = dict()
with open('../genome/disease_ppi.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease_ppi[str2[0]] = set(str2[1].split(';'))
    infile.close()

disease_pathway = dict()
with open('../genome/disease_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease_pathway[str2[0]] = set(str2[1].split(';'))
    infile.close()

comorbidity_snp = dict()
with open('../overlap/comorbidity_snp.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        comorbidity_snp[(str2[0], str2[1])] = set(str2[4].split(';'))
    infile.close()


comorbidity_gene = dict()
with open('../overlap/comorbidity_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        comorbidity_gene[(str2[0], str2[1])] = set(str2[4].split(';'))
    infile.close()


comorbidity_ppi = dict()
with open('../overlap/comorbidity_ppi.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        comorbidity_ppi[(str2[0], str2[1])] = set(str2[4].split(';'))
    infile.close()

comorbidity_pathway = dict()
with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        comorbidity_pathway[(str2[0], str2[1])] = set(str2[4].split(';'))
    infile.close()

list1 = list()
for each in comorbidity_snp:
    code1, code2 = each
    a = len(comorbidity_snp[each])
    b = len(disease_snp[code1]) - a
    c = len(disease_snp[code2]) - a
    d = all_variant_count - a - b - c
    [odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
    list1.append([code1, code2, odds, p])
df = pd.DataFrame(list1, columns=['code1', 'code2', 'odds', 'p'])
df['q'] = multitest.fdrcorrection(df['p'])[1]
df = df[df['q']<0.05]
print('snp significance: ' + str(df.shape[0]))
snp_significant = set()
for each in df.values.tolist():
    snp_significant.add((each[0], each[1]))


list1 = list()
for each in comorbidity_gene:
    code1, code2 = each
    a = len(comorbidity_gene[each])
    b = len(disease_gene[code1]) - a
    c = len(disease_gene[code2]) - a
    d = all_gene_count - a - b - c
    [odds, p] = fisher_exact([[a,b], [c,d]], alternative='greater')
    list1.append([code1, code2, odds, p])
df = pd.DataFrame(list1, columns=['code1', 'code2', 'odds', 'p'])
df['q'] = multitest.fdrcorrection(df['p'])[1]
df = df[df['q']<0.05]
print('gene significance: ' + str(df.shape[0]))
gene_significant = set()
for each in df.values.tolist():
    gene_significant.add((each[0], each[1]))


result = dict()
for each in comorbidity_ppi:
    result[each] = 0
i = 0
while i < 10000:
    print(i)
    disease_gene1 = dict()
    disease_ppi1 = dict()
    for each in disease_gene:
        set1 = set(random.sample(all_gene, len(disease_gene[each])))
        disease_gene1[each] = set1
    for each in disease_gene1:
        set1 = disease_gene1[each]
        set3 = set()
        for each1 in set1:
            if each1 in ppi_dict:
                set2 = ppi_dict[each1]
                for each2 in set2:
                    if each1 < each2:
                        set3.add(each1 + '~' + each2)
                    else:
                        set3.add(each2 + '~' + each1)
        disease_ppi1[each] = set3
    for each in comorbidity_ppi:
        a = len(comorbidity_ppi[each])
        code1, code2 = each
        set1 = set()
        if (code1 in disease_ppi1) & (code2 in disease_ppi1):
            set2 = disease_ppi1[code1] & disease_ppi1[code2]
            for each1 in set2:
                gene1, gene2 = each1.split('~')
                if len(set([gene1, gene2]) & disease_gene1[code1]) == 2:
                    continue
                if len(set([gene1, gene2]) & disease_gene1[code2]) == 2:
                    continue
                if (gene1 in disease_gene1[code1]) & (gene2 in disease_gene1[code2]):
                    set1.add(each1)
                elif (gene2 in disease_gene1[code1]) & (gene1 in disease_gene1[code2]):
                    set1.add(each1)
        if len(set1) >= a:
            result[each] += 1
    i += 1

list1 = list()
for each in result:
    list1.append([each[0], each[1], result[each]/10000])
df = pd.DataFrame(list1, columns=['code1', 'code2', 'p'])
df['q'] = multitest.fdrcorrection(df['p'])[1]
df = df[df['q']<0.05]
print('ppi significance: ' + str(df.shape[0]))
ppi_significant = set()
for each in df.values.tolist():
    ppi_significant.add((each[0], each[1]))


result = dict()
for each in comorbidity_pathway:
    result[each] = 0
i = 0
while i < 10000:
    print(i)
    disease_gene1 = dict()
    disease_pathway1 = dict()
    for each in disease_gene:
        set1 = set(random.sample(all_gene, len(disease_gene[each])))
        disease_gene1[each] = set1
    list2, pathway_gene = pathwayEnrich(disease_gene1)
    df = pd.DataFrame(list2, columns=['disease', 'pathway', 'odds', 'p'])
    df = df[df['odds']>1]
    df['q'] = multitest.fdrcorrection(df['p'])[1]
    df = df[df['q']<0.05]
    list3 = df.values.tolist()
    for each in list3:
        if each[0] not in disease_pathway1:
            disease_pathway1[each[0]] = set()
        disease_pathway1[each[0]].add(each[1])
    for each in comorbidity_pathway:
        a = len(comorbidity_pathway[each])
        code1, code2 = each
        set1 = set()
        if code1 in disease_pathway1:
            set2 = disease_pathway1[code1]
            for each1 in set2:
                set3 = pathway_gene[each1]
                if len(disease_gene1[code2] & set3) != 0:
                    set1.add(each1)
        if code2 in disease_pathway1:
            set2 = disease_pathway1[code2]
            for each1 in set2:
                set3 = pathway_gene[each1]
                if len(disease_gene1[code1] & set3) != 0:
                    set1.add(each1)
        if len(set1) >= a:
            result[each] += 1
    i += 1

list1 = list()
for each in result:
    list1.append([each[0], each[1], result[each]])
df = pd.DataFrame(list1, columns=['code1', 'code2', 'p'])
df['q'] = multitest.fdrcorrection(df['p'])[1]
df = df[df['q']<0.05]
print('pathway significance: ' + str(df.shape[0]))

pathway_significant = set()
for each in df.values.tolist():
    pathway_significant.add((each[0], each[1]))


all_interpreted_comorbidity = set()
for each in comorbidity_snp:
    all_interpreted_comorbidity.add((each[0], each[1]))
for each in comorbidity_gene:
    all_interpreted_comorbidity.add((each[0], each[1]))
for each in comorbidity_ppi:
    all_interpreted_comorbidity.add((each[0], each[1]))
for each in comorbidity_pathway:
    all_interpreted_comorbidity.add((each[0], each[1]))

list1 = list()
for each in all_interpreted_comorbidity:
    if each in snp_significant:
        snp = 'Y'
    else:
        snp = ''
    if each in gene_significant:
        gene = 'Y'
    else:
        gene = ''
    if each in ppi_significant:
        ppi = 'Y'
    else:
        ppi = ''
    if each in pathway_significant:
        pathway = 'Y'
    else:
        pathway = ''
    list1.append([each[0], each[1], snp, gene, ppi, pathway])
df = pd.DataFrame(list1, columns=['Code1', 'Code2', 'Significance of SNP', 'Significance of gene',
                                  'Significance of PPI', 'Significance of pathway'])
df.to_csv('a.csv', index=False)