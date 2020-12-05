import networkx as nx
from networkx.algorithms import sigma


if __name__ == '__main__':

    set1 = set()
    print('snp and gene')
    with open('../overlap/comorbidity_snp.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            if (code1 == 'C64') & (code2 == 'I95'):
                continue
            set1.add((code1, code2))
        infile.close()

    with open('../overlap/comorbidity_gene.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            if (code1 == 'C64') & (code2 == 'I95'):
                continue
            set1.add((code1, code2))
        infile.close()

    G = nx.Graph()
    G.add_edges_from(set1)
    sig_ma = sigma(G, niter=100, nrand=10)
    print(sig_ma)



    set1 = set()
    print('ppi and pahtway')
    with open('../overlap/comorbidity_ppi.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            set1.add((code1, code2))
        infile.close()

    with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1 = str2[0]
            code2 = str2[1]
            set1.add((code1, code2))
        infile.close()

    G = nx.Graph()
    G.add_edges_from(set1)
    sig_ma = sigma(G, niter=100, nrand=10)
    print(sig_ma)