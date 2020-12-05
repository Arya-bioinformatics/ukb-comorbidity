import community
import  networkx as nx
from community import community_louvain
import pandas as pd

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

disease_class = dict()
disease_name = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    code = each[0]
    name = each[1]
    c = each[2]
    disease_class[code] = c
    disease_name[code] = name

# ----------------------------- module detection (LG-Network) --------------------------- #
G = nx.Graph()
for each in comorbidity_loci:
    code1, code2 = each
    if code1 in top_025_disease_loci:
        continue
    if code2 in top_025_disease_loci:
        continue
    G.add_edge(code1, code2, weight=1.0)

partition = community_louvain.best_partition(G, resolution=1, weight='weight')
print(len(set(partition.values())))

mod = community.modularity(partition, G, weight='weight')
print(mod)

list1 = list()
for each in partition:
    name = disease_name[each]
    c = disease_class[each]
    list1.append([each, name, c, partition[each]])
df = pd.DataFrame(list1, columns=['node', 'name', 'category', 'cluster'])
df.to_csv('LG_network_module.csv', index=False)


# ----------------------------- module detection (NG-Network) --------------------------- #
G = nx.Graph()
for each in comorbidity_network:
    code1, code2 = each
    if code1 in top_025_disease_network:
        continue
    if code2 in top_025_disease_network:
        continue
    G.add_edge(code1, code2, weight=1.0)

partition = community_louvain.best_partition(G, resolution=1, weight='weight')

mod = community.modularity(partition, G, weight='weight')
print(mod)

list1 = list()
dict1 = dict()
for each in partition:
    name = disease_name[each]
    c = disease_class[each]
    list1.append([each, name, c, partition[each]])
    if partition[each] not in dict1:
        dict1[partition[each]] = set()
    dict1[partition[each]].add(each)
df = pd.DataFrame(list1, columns=['node', 'name', 'category', 'cluster'])
df.to_csv('NG_network_module.csv', index=False)
print('total modules: ' + str(len(dict1)))
for each in dict1:
    print([each, len(dict1[each])])