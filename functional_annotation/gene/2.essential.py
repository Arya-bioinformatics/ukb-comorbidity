from scipy.stats import fisher_exact
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


all_gene= set()
with open('../genome/grch37_download/genome_assemblies/gene_assembly_grch37.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\n')
        str2 = str1.split('\t')
        gene_id = str2[0]
        if type == 'pseudo':
            continue
        all_gene.add(gene_id)
    infile.close()


lethality_phenotype = set()
with open('../overlap/MGI_download/VOC_MammalianPhenotype.rpt.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        if 'lethal' in line:
            if ('embryonic' in line) | ('prenatal' in line) | ('postnatal' in line):
                lethality_phenotype.add(str2[0])
    infile.close()

gene_lethality = set()
with open('../overlap/MGI_download/MGI_PhenoGenoMP.rpt.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        if not str2[3].startswith('MP'):
            continue
        phenotype = str2[3]
        if phenotype not in lethality_phenotype:
            continue
        gene = str2[5]
        gene_lethality.add(gene)
    infile.close()

gene_valid_knockout = set()
with open('../overlap/MGI_download/MGI_PhenotypicAllele.rpt.txt', 'r') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        gene = str2[6]
        if gene not in gene_lethality:
            continue
        type = str2[4]
        if 'knockout' not in type:
            continue
        phenotype = str2[10]
        if not phenotype.startswith('MP'):
            continue
        set1 = set(phenotype.split(','))
        # 'MP:0010768' - mortality/aging	the observable characteristics related to the ability of a mammalian organism to live and age that are manifested throughout development and life span
        # this is due to hierarchical structure
        if 'MP:0010768' not in set1:
            continue
        gene_valid_knockout.add(gene)
    infile.close()

print(len(gene_valid_knockout))

essential_gene = set()
with open('../overlap/MGI_download/HMD_HumanPhenotype.rpt.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        entrezid = str2[1].strip(' ')
        mgi_id = str2[5].strip(' ')
        if mgi_id not in gene_valid_knockout:
            continue
        phenotype = str2[6].strip(' ')
        set1 = set(phenotype.split(' '))
        if 'MP:0010768' not in set1:
            continue
        essential_gene.add(entrezid)
    infile.close()

print('total essential gene:' + str(len(essential_gene)))

for flag in ['with', 'without']:
    print('\n' + flag + '\n')
    if flag == 'with':
        path1 = '../genome/disease_gene.txt'
        path2 = '../overlap/comorbidity_gene.txt'
        path3 = 'a.pdf'
        path4 = 'b.pdf'
        # continue
    if flag == 'without':
        path1 = '../genome/disease_gene_rmhla.txt'
        path2 = '../overlap/comorbidity_gene_rmhla.txt'
        path3 = 'a1.pdf'
        path4 = 'b1.pdf'
        # continue

    disease_gene = set()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            disease_gene |= set(str2[1].split(';'))
        infile.close()

    comorbidity_gene = set()
    with open(path2, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            comorbidity_gene |= set(str2[4].split(';'))
        infile.close()

    print(' --------- gene count ---------')
    print('all gene count: ' + str(len(all_gene)))
    print('disease gene count: ' + str(len(disease_gene)))
    print('comorbidity gene count: ' + str(len(comorbidity_gene)))
    print('essential gene count: ' + str(len(essential_gene)))


    print('# ----- essential gene overlap ----- #')
    print(len(comorbidity_gene & essential_gene) / float(len(comorbidity_gene)))
    print(len(disease_gene & essential_gene) / float(len(disease_gene)))

    print(' --------- comorbidity enrich ---------- ')
    a = len(comorbidity_gene & essential_gene)
    b = len(essential_gene) - a
    c = len(comorbidity_gene) - a
    d = len(all_gene) - a - b - c
    [odds, p] = fisher_exact([[a, b], [c, d]])
    print([odds, p])

    print('--------- comorbidity enrich (backgroud disease) ----------')
    a = len(comorbidity_gene & essential_gene)
    b = len(essential_gene & disease_gene) - a
    c = len(comorbidity_gene) - a
    d = len(disease_gene) - a - b - c
    [odds, p] = fisher_exact([[a, b], [c, d]])
    print([odds, p])

    print('--------- disease enrich ----------')
    a = len(disease_gene & essential_gene)
    b = len(essential_gene) - a
    c = len(disease_gene) - a
    d = len(all_gene) - a - b - c
    [odds, p] = fisher_exact([[a, b], [c, d]])
    print([odds, p])

    fig = plt.figure(figsize=(10, 6))
    c = venn2([comorbidity_gene, essential_gene], ('Comorbidity-genes', 'Essential genes'), alpha=1)
    c.get_patch_by_id('10').set_color('#F8766D')
    c.get_patch_by_id('01').set_color('#00BA38')
    for text in c.set_labels:
        text.set_fontsize(30)
    for text in c.subset_labels:
        text.set_fontsize(30)
    c.get_patch_by_id('10').set_edgecolor('black')
    c.get_patch_by_id('01').set_edgecolor('black')
    plt.savefig(path3, bbox_inches="tight")
    plt.close()

    fig = plt.figure(figsize=(10, 6))
    c = venn2([disease_gene, essential_gene], ('Disease-genes', 'Essential genes'), alpha=1)
    c.get_patch_by_id('10').set_color('#F8766D')
    c.get_patch_by_id('01').set_color('#00BA38')
    for text in c.set_labels:
        text.set_fontsize(30)
    for text in c.subset_labels:
        text.set_fontsize(30)
    c.get_patch_by_id('10').set_edgecolor('black')
    c.get_patch_by_id('01').set_edgecolor('black')
    plt.savefig(path4, bbox_inches="tight")
    plt.close()