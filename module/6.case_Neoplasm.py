import pandas as pd


disease_class = dict()
disease_name = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    disease_class[each[0]] = each[2]
    disease_name[each[0]] = each[1]

cluster_disease = dict()
df = pd.read_csv('../module/LG_network_module_Q0.4007486988.csv')
for each in df.values.tolist():
    disease, cluster = each[0], each[3]
    if cluster not in cluster_disease:
        cluster_disease[cluster] = set()
    cluster_disease[cluster].add(disease)

for each in cluster_disease:
    if 'C44' in cluster_disease[each]:
        der_neop_cluster_disease = cluster_disease[each]
    if 'N40' in cluster_disease[each]:
        male_neop_urin_cluster_disease = cluster_disease[each]


comorbidity_gene = dict()
df = pd.read_table('../overlap/comorbidity_gene.txt')
for each in df.values.tolist():
    code1, code2 = each[:2]
    c1 = disease_class[code1]
    c2 = disease_class[code2]
    if 'Neoplasm' in set([c1, c2]):
        comorbidity_gene[(code1, code2)] = set(each[5].split(';'))

list1 = list()
list2 = list()
cytoscape_input = dict()
for each in comorbidity_gene:
    if len(set(each) & der_neop_cluster_disease) == 2:
        list1 += comorbidity_gene[each]
        cytoscape_input[each] = comorbidity_gene[each]
    if len(set(each) & male_neop_urin_cluster_disease) == 2:
        list2 += comorbidity_gene[each]
        cytoscape_input[each] = comorbidity_gene[each]

der_neop_gene_count = list()
for each in set(list1):
    der_neop_gene_count.append([each, list1.count(each)])

male_neop_urin_gene_count = list()
for each in set(list2):
    male_neop_urin_gene_count.append([each, list2.count(each)])

df1 = pd.DataFrame(der_neop_gene_count, columns=['gene', 'count'])
df2 = pd.DataFrame(male_neop_urin_gene_count, columns=['gene', 'count'])
df1 = df1.sort_values(by='count', ascending=False)
df2 = df2.sort_values(by='count', ascending=False)

print(len(set(df1['gene'])))
print(len(set(df2['gene'])))
print(set(df1['gene']) & set(df2['gene']))

edges = set()
nodes = set()
for each in cytoscape_input:
    code1, code2 = each
    edges.add((code1, code2))
    nodes.add((code1, 'disease', disease_name[code1], disease_class[code1]))
    nodes.add((code2, 'disease', disease_name[code2], disease_class[code2]))
    set1 = cytoscape_input[each]
    for each1 in set1:
        edges.add((code1, each1))
        edges.add((code2, each1))
        nodes.add((each1, 'gene', '', ''))

with open('edges1.txt', 'w+') as outfile:
    outfile.write('node1\tnode2\n')
    for each in edges:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()

with open('node_atribute1.txt', 'w+') as outfile:
    outfile.write('node\ttype\tname\tcategory\n')
    for each in nodes:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()