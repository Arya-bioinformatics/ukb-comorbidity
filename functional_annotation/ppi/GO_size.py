
list1 = list()
with open('../overlap/GO_geneset/download/c5.bp.v7.0.entrez.gmt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        bp = str2[0]
        set1 = set(str2[2:])
        if '' in set1:
            set1.remove('')
        if len(set1) == 0:
            continue
        list1.append(len(set1))
    infile.close()

list1.sort()

print(len(list1))
print(len([i for i in list1 if i < 200]))

# plt.plot(range(0, len(list1)), list1)
# plt.show()