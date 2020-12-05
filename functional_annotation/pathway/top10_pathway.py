import pandas as pd
from matplotlib import pyplot as plt
from collections import Counter


for flag in ['with', 'without']:
    print('\n' + flag + '\n')
    if flag == 'with':
        path1 = '../overlap/comorbidity_pathway.txt'
        path2 = 'a.pdf'
    if flag == 'without':
        path1 = '../overlap/comorbidity_pathway_rmhla.txt'
        path2 = 'a1.pdf'

    comorbidity_pathway = dict()
    pathway_comorbidity_count = dict()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            set1 = set(str2[4].split(';'))
            for each in set1:
                if each not in pathway_comorbidity_count:
                    pathway_comorbidity_count[each] = 0
                pathway_comorbidity_count[each] += 1
            comorbidity_pathway[(code1, code2)] = set1
        infile.close()


    list1 = list()
    for each in pathway_comorbidity_count:
        list1.append([each, pathway_comorbidity_count[each]])

    df = pd.DataFrame(list1, columns=['pathway', 'count'])
    df = df.sort_values(by='count', ascending=False)

    df1 = df[0:10]
    list1 = df1.values.tolist()
    list2 = list()
    list3 = list()
    for each in list1:
        pathway, count = each
        list2.append(pathway)
        list3.append(count)

    for each in list2:
        print(each)

    plt.figure(figsize=(6, 4))
    plt.bar(range(0, len(list3)), list3, alpha=0.5, color='steelblue', edgecolor='black', linewidth=0.5, width=0.9)
    plt.xticks(range(0, len(list2)), list2, rotation=90)
    plt.xlabel('Pathway', fontsize=15)
    plt.ylabel('Comorbidity count', fontsize=15)
    plt.savefig(path2, bbox_inches='tight')

    set3 = set()
    for each in comorbidity_pathway:
        set2 = comorbidity_pathway[each]
        if len(set(list2) & set2) != 0:
            set3.add(each)
    list1 = list()
    for each in set3:
        list1.append(each[0])
        list1.append(each[1])
    x = Counter(list1)


    print('------------- total comorbidity with pathway overlap -----------')
    print(len(comorbidity_pathway))

    print('--------- total comorbidty with top 10 pathway overlap ---------')
    print(len(set3))

    print('---------------------------ratio  ------------------------------')
    print(len(set3)/len(comorbidity_pathway))