import pandas as pd


disease_class = dict()
disease_name = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    disease_name[each[0]] = each[1]
    disease_class[each[0]] = each[2]


# ----------------------------- find top disease --------------------------- #
top = 0.25

comorbidity_loci = set()
with open('../overlap/comorbidity_snp.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_loci.add((code1, code2))
    infile.close()

with open('../overlap/comorbidity_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_loci.add((code1, code2))
    infile.close()

disease_degree = dict()
for each in comorbidity_loci:
    code1, code2 = each
    if code1 not in disease_degree:
        disease_degree[code1] = 0
    disease_degree[code1] += 1
    if code2 not in disease_degree:
        disease_degree[code2] = 0
    disease_degree[code2] += 1

list1 = list()
for each in disease_degree:
    list1.append([each, disease_degree[each]])
df = pd.DataFrame(list1, columns=['disease', 'degree'])
top_025_disease_loci = set(df.loc[df['degree']>df.shape[0]*top, 'disease'])
print('total top 25% disease loci: ' + str(len(top_025_disease_loci)))

comorbidity_network = set()
with open('../overlap/comorbidity_ppi.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_network.add((code1, code2))
    infile.close()

with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        comorbidity_network.add((code1, code2))
    infile.close()

disease_degree = dict()
for each in comorbidity_network:
    code1, code2 = each
    if code1 not in disease_degree:
        disease_degree[code1] = 0
    disease_degree[code1] += 1
    if code2 not in disease_degree:
        disease_degree[code2] = 0
    disease_degree[code2] += 1

list1 = list()
for each in disease_degree:
    list1.append([each, disease_degree[each]])
df = pd.DataFrame(list1, columns=['disease', 'degree'])
top_025_disease_network = set(df.loc[df['degree'] > df.shape[0] * top, 'disease'])
print('total top 25% disease network: ' + str(len(top_025_disease_network)))

disease_module_loci = dict()
df = pd.read_csv('../module/LG_network_module_Q0.4007486988.csv')
list1 = df.values.tolist()
for each in list1:
    disease_module_loci[each[0]] = each[3]
for each in top_025_disease_loci:
    disease_module_loci[each] = -1

disease_module_network = dict()
df = pd.read_csv('../module/NG_network_module_Q0.3200458.csv')
list1 = df.values.tolist()
for each in list1:
    disease_module_network[each[0]] = each[3]
for each in top_025_disease_network:
    disease_module_network[each] = -1

edges = set()
nodes = set()
for each in comorbidity_loci:
    code1, code2 = each
    if (code1 in disease_module_loci) & (code2 in disease_module_loci):
        c1 = disease_module_loci[code1]
        c2 = disease_module_loci[code2]
        if c1 == c2:
            edges.add((code1, code2, 'inner'))
        else:
            edges.add((code1, code2, 'outer'))
        nodes.add((code1, disease_name[code1], disease_class[code1], str(disease_module_loci[code1])))
        nodes.add((code2, disease_name[code2], disease_class[code2], str(disease_module_loci[code2])))

with open('edges4.txt', 'w+') as outfile:
    outfile.write('node1\tnode2\tcluster\n')
    for each in edges:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()

with open('node_atribute4.txt', 'w+') as outfile:
    outfile.write('node\tdescription\tcategory\tmodule\n')
    for each in nodes:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()

edges = set()
nodes = set()
for each in comorbidity_network:
    code1, code2 = each
    if (code1 in disease_module_network) & (code2 in disease_module_network):
        c1 = disease_module_network[code1]
        c2 = disease_module_network[code2]
        if c1 == c2:
            edges.add((code1, code2, 'inner'))
        else:
            edges.add((code1, code2, 'outer'))
        nodes.add((code1, disease_name[code1], disease_class[code1], str(disease_module_network[code1])))
        nodes.add((code2, disease_name[code2], disease_class[code2], str(disease_module_network[code2])))

with open('edges5.txt', 'w+') as outfile:
    outfile.write('node1\tnode2\tcluster\n')
    for each in edges:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()

with open('node_atribute5.txt', 'w+') as outfile:
    outfile.write('node\tdescription\tcategory\tmodule\n')
    for each in nodes:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()