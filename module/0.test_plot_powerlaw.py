import networkx as nx
from matplotlib import pyplot as plt
import pandas as pd
import math
import powerlaw
import statsmodels.api as sm
from collections import Counter


if __name__ == '__main__':

    comorbidity_snp = set()
    with open('../overlap/comorbidity_snp.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            comorbidity_snp.add((code1, code2))
        infile.close()

    comorbidity_gene = set()
    with open('../overlap/comorbidity_gene.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            comorbidity_gene.add((code1, code2))
        infile.close()

    comorbidity_ppi = set()
    with open('../overlap/comorbidity_ppi.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            comorbidity_ppi.add((code1, code2))
        infile.close()

    comorbidity_pathway = set()
    with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            comorbidity_pathway.add((code1, code2))
        infile.close()

    print('# ---------------- test and plot snp-gene degree distribution ----------------- #')
    G = nx.Graph()
    for each in (comorbidity_snp | comorbidity_gene):
        code1, code2 = each
        G.add_edge(code1, code2, weight=1.0)

    df = pd.DataFrame(G.degree(), columns=['disease', 'degree'])
    df = df.sort_values(by='degree', ascending=True)
    df.to_csv('a.csv', index=False)
    list1 = list(df['degree'])

    # -------- 检验是否服从幂律分布 alfa在2-3之间 ------------ #
    results = powerlaw.Fit(list1)
    print(results.power_law.alpha)
    print(results.power_law.xmin)
    Xmin = results.power_law.xmin
    list2 = list()
    for each in list1:
        if each < Xmin:
            continue
        list2.append(each)
    degree_count = Counter(list2)
    df = pd.DataFrame(degree_count.items(), columns=['degree', 'cnt'])
    l1 = list(df['degree'])
    l2 = list(df['cnt'])
    X = [math.log(i, 10) for i in l1]
    y = [math.log(i, 10) for i in l2]
    est = sm.OLS(y, X)
    est = est.fit()
    print(est.summary())


    dict1 = dict()
    for each in list1:
        if each not in dict1:
            dict1[each] = 0
        dict1[each] += 1

    list1 = list()
    for each in dict1:
        count = dict1[each]
        list1.append([math.log(each), math.log(count)])

    df = pd.DataFrame(list1, columns=['degree', 'count'])
    df = df.sort_values(by='degree', ascending=True)
    plt.figure(figsize=(6,4))
    plt.subplot(121)
    plt.scatter(df['degree'], df['count'])
    plt.xlabel('Log10(degree)')
    plt.ylabel('Log10(num)')
    plt.title('Loci level')




    print('# ---------------- test and plot ppi-pathway degree distribution ----------------- #')
    G = nx.Graph()
    for each in (comorbidity_ppi | comorbidity_pathway):
        code1, code2 = each
        G.add_edge(code1, code2, weight=1.0)

    df = pd.DataFrame(G.degree(), columns=['disease', 'degree'])
    df = df.sort_values(by='degree', ascending=True)
    df.to_csv('a.csv', index=False)
    list1 = list(df['degree'])

    # -------- 检验是否服从幂律分布 alfa在2-3之间 ------------ #
    results = powerlaw.Fit(list1)
    print(results.power_law.alpha)
    print(results.power_law.xmin)
    Xmin = results.power_law.xmin
    list2 = list()
    for each in list1:
        if each < Xmin:
            continue
        list2.append(each)
    degree_count = Counter(list2)
    df = pd.DataFrame(degree_count.items(), columns=['degree', 'cnt'])
    l1 = list(df['degree'])
    l2 = list(df['cnt'])
    X = [math.log(i, 10) for i in l1]
    y = [math.log(i, 10) for i in l2]
    est = sm.OLS(y, X)
    est = est.fit()
    print(est.summary())

    dict1 = dict()
    for each in list1:
        if each not in dict1:
            dict1[each] = 0
        dict1[each] += 1

    list1 = list()
    for each in dict1:
        count = dict1[each]
        list1.append([math.log(each), math.log(count)])

    df = pd.DataFrame(list1, columns=['degree', 'count'])
    df = df.sort_values(by='degree', ascending=True)
    plt.subplot(122)
    plt.scatter(df['degree'], df['count'])
    plt.xlabel('Log10(degree)')
    plt.ylabel('Log10(num)')
    plt.title('Network level')
    # plt.show()
    plt.savefig('a.pdf', bbox_inches='tight')