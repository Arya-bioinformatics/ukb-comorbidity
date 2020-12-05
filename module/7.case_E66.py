import pandas as pd
from itertools import combinations
from collections import Counter


disease_class = dict()
disease_name = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    disease_class[each[0]] = each[2]
    disease_name[each[0]] = each[1]

cluster_disease = dict()
df = pd.read_csv('../module/NG_network_module_Q0.3200458.csv')
list1 = df.values.tolist()
for each in list1:
    disease, name, category, cluster = each
    if cluster not in cluster_disease:
        cluster_disease[cluster] = set()
    cluster_disease[cluster].add(disease)
    if disease == 'M07':
        joint_mascular_neur_spine_cluster = cluster
    if disease == 'K80':
        hepatobiliary_cluster = cluster

comorbidity_pathway = dict()
df = pd.read_table('../overlap/comorbidity_pathway.txt')
list1 = df.values.tolist()
for each in list1:
    comorbidity_pathway[(each[0], each[1])] = set(each[4].split(';'))

print('# ------------------ Joint-Muscular-Neurological-Spine -------------- #')
list1 = list()
for each in comorbidity_pathway:
    if ('E66' in each) & (len(set(each) & cluster_disease[joint_mascular_neur_spine_cluster]) !=0):
        list1 += list(comorbidity_pathway[each])
for each in list1:
    print(each)

print('# ------------------------ Hepatobiliary pancreas --------------------- #')
list1 = list()
for each in comorbidity_pathway:
    if ('E66' in each) & (len(set(each) & cluster_disease[hepatobiliary_cluster]) !=0):
        list1 += list(comorbidity_pathway[each])
for each in list1:
    print(each)

cytoscape_input = dict()
list1 = list()
dict1 = dict()
for each in comorbidity_pathway:
    if ('E66' in each) & (len(set(each) & cluster_disease[joint_mascular_neur_spine_cluster]) !=0):
        list1 += list(comorbidity_pathway[each])
        dict1[each] = comorbidity_pathway[each]
pathway_count_1 = Counter(list1)
# in order to include all the pathway with same frequency, set 6
top5_pathway_1 = set([i[0] for i in pathway_count_1.most_common(6)])

list2 = list()
dict2 = dict()
for each in comorbidity_pathway:
    if ('E66' in each) & (len(set(each) & cluster_disease[hepatobiliary_cluster]) !=0):
        list2 += list(comorbidity_pathway[each])
        dict2[each] = comorbidity_pathway[each]
pathway_count_2 = Counter(list2)
# in order to include all the pathway with same frequency, set 6
top5_pathway_2 = set([i[0] for i in pathway_count_2.most_common(10)])

edges = set()
nodes = set()
for each in dict1:
    code1, code2 = each
    set1 = dict1[each] & top5_pathway_1
    if len(set1) != 0:
        edges.add((code1, code2))
        nodes.add((code1, 'disease', disease_name[code1], disease_class[code1]))
        nodes.add((code2, 'disease', disease_name[code2], disease_class[code2]))
    for each1 in set1:
        edges.add((code1, each1))
        edges.add((code2, each1))
        nodes.add((each1, 'pathway', '', ''))
for each in dict2:
    code1, code2 = each
    set1 = dict2[each] & top5_pathway_2
    if len(set1) != 0:
        edges.add((code1, code2))
        nodes.add((code1, 'disease', disease_name[code1], disease_class[code1]))
        nodes.add((code2, 'disease', disease_name[code2], disease_class[code2]))
    for each1 in set1:
        edges.add((code1, each1))
        edges.add((code2, each1))
        nodes.add((each1, 'pathway', '', ''))

with open('edges3.txt', 'w+') as outfile:
    outfile.write('node1\tnode2\n')
    for each in edges:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()

with open('node_atribute3.txt', 'w+') as outfile:
    outfile.write('node\ttype\tname\tcategory\n')
    for each in nodes:
        outfile.write('\t'.join(each) + '\n')
    outfile.close()