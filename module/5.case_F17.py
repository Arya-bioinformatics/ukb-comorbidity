import pandas as pd
from itertools import combinations

cluster_disease = dict()
disease_name = dict()
disease_class = dict()
df = pd.read_csv('../module/LG_network_module_Q0.4007486988.csv')
list1 = df.values.tolist()
for each in list1:
    disease, name, category, cluster = each
    disease_name[disease] = name
    disease_class[disease] = category
    if cluster not in cluster_disease:
        cluster_disease[cluster] = set()
    cluster_disease[cluster].add(disease)
    if disease == 'F17':
        F17_cluster = cluster

comorbidity_gene = dict()
df = pd.read_table('../overlap/comorbidity_gene.txt')
list1 = df.values.tolist()
for each in list1:
    comorbidity_gene[(each[0], each[1])] = set(each[5].split(';'))

all_comb = combinations(cluster_disease[F17_cluster], 2)
diseasePair_set = set()
for each in all_comb:
    code1, code2 = each
    if code1 < code2:
        diseasePair_set.add((code1, code2))
    else:
        diseasePair_set.add((code2, code1))
module_comorbidity = set(comorbidity_gene.keys()) & diseasePair_set

list1 = list()
for each in module_comorbidity:
    if 'F17' in each:
        print([disease_name[each[0]], disease_name[each[1]], disease_class[each[0]], disease_class[each[1]]])

dict1 = dict()
for each in comorbidity_gene:
    if ('F17' in each) & (each in module_comorbidity):
        for each1 in comorbidity_gene[each]:
            if each1 not in dict1:
                dict1[each1] = 0
            dict1[each1] += 1

list1 = list()
for each in dict1:
    list1.append([each, dict1[each]])

df = pd.DataFrame(list1, columns=['gene', 'count'])
df = df.sort_values(by='count', ascending=False)


cytoscape_input = dict()
for each in comorbidity_gene:
    if ('F17' in each) & (len(set(each) & cluster_disease[F17_cluster]) == 2):
        cytoscape_input[each] = comorbidity_gene[each]

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

with open('edges2.txt', 'w+') as outfile:
    outfile.write('node1\tnode2\n')
    for each in edges:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()

with open('node_atribute2.txt', 'w+') as outfile:
    outfile.write('node\ttype\tname\tcategory\n')
    for each in nodes:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()
