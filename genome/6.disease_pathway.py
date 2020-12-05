from scipy.stats import fisher_exact
from statsmodels.stats import multitest
import pandas as pd


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
    return list1


for flag in ['with', 'without']:
    if flag == 'with':
        path1 = '../genome/disease_gene.txt'
        path2 = '../genome/disease_pathway.txt'
    elif flag == 'without':
        path1 = '../genome/disease_gene_rmhla.txt'
        path2 = '../genome/disease_pathway_rmhla.txt'

    disease_gene = dict()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_gene[str2[0]] = set(str2[1].split(';'))
        infile.close()

    result = pathwayEnrich(disease_gene)
    df = pd.DataFrame(result, columns=['disease', 'pathway', 'odds', 'p'])
    df = df[df['odds']>1]
    df['q'] = multitest.fdrcorrection(df['p'])[1]
    df = df[df['q']<0.05]
    all_disease = set(df['disease'])
    disease_pathway = dict()
    for each in all_disease:
        disease_pathway[each] = set(df.loc[df['disease']==each, 'pathway'])

    with open(path2, 'w+') as outfile:
        outfile.write('disease\tpathway\n')
        for each in disease_pathway:
            outfile.write(each + '\t' + ';'.join(disease_pathway[each]) + '\n')
        outfile.close()