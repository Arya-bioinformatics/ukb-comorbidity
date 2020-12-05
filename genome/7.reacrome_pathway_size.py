import matplotlib.pyplot as plt
import numpy as np

kegg = list()
biocarta = list()
reactome = list()
pid = list()

pathway_size = dict()
with open('../genome/MSigDB/c2.cp.v6.2.entrez.gmt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        if str1.startswith('KEGG'):
            set1 = set(str2[2:])
            kegg.append(len(set1))
        elif str1.startswith('BIOCARTA'):
            set1 = set(str2[2:])
            biocarta.append(len(set1))
        elif str1.startswith('REACTOME'):
            set1 = set(str2[2:])
            reactome.append(len(set1))
        elif str1.startswith('PID'):
            set1 = set(str2[2:])
            pid.append(len(set1))
        pathway_size[str.lower(str2[0])] = len(str2[2:])
    infile.close()

kegg.sort()
biocarta.sort()
reactome.sort()
pid.sort()
print('total number of pathway: reactome, kegg, biocarta, pid')
print(len(reactome))
print(len(kegg))
print(len(biocarta))
print(len(pid))
print('smallest size: reactome, kegg, biocarta, pid')
print(min(reactome))
print(min(kegg))
print(min(biocarta))
print(min(pid))
print('maximum size: reactome, kegg, biocarta, pid')
print(max(reactome))
print(max(kegg))
print(max(biocarta))
print(max(pid))
print('avarage size: reactome, kegg, biocarta, pid')
print(np.mean(reactome))
print(np.mean(kegg))
print(np.mean(biocarta))
print(np.mean(pid))


plt.figure(figsize=(16, 3))
plt.subplot(1,4,1)
plt.hist(reactome, color='blue')
plt.xlabel('Patyway size\nmin:9 max:933 mean:56')
plt.ylabel('Number of pathway')
plt.title('Reactome')

plt.subplot(1,4,2)
plt.hist(kegg, color='blue')
plt.xlabel('Patyway size\nmin:10 max:389 mean:69')
plt.ylabel('Number of pathway')
plt.title('KEGG')

plt.subplot(1,4,3)
plt.hist(biocarta, color='blue')
plt.xlabel('Patyway size\nmin:6 max:87 mean:21')
plt.ylabel('Number of pathway')
plt.title('BioCarta')

plt.subplot(1,4,4)
plt.hist(pid, color='blue')
plt.xlabel('Patyway size\nmin:10 max:137 mean:41')
plt.ylabel('Number of pathway')
plt.title('PID')
# plt.show()
plt.savefig('a.pdf', bbox_inches='tight')

# after removing large pathway
disease_pathway = set()
with open('../genome/disease_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease_pathway |= set(str2[1].split(';'))
    infile.close()

comorbidity_pathway = set()
with open('../overlap/comorbidity_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        comorbidity_pathway |= set(str2[4].split(';'))
    infile.close()

list1 = list()
list2 = list()
list3 = list()
list4 = list()
for each in disease_pathway:
    size = pathway_size[each]
    if each.startswith('reactome'):
        list1.append(size)
    else:
        list2.append(size)
for each in comorbidity_pathway:
    size = pathway_size[each]
    if each.startswith('reactome'):
        list3.append(size)
    else:
        list4.append(size)

list1.sort()
list2.sort()
list3.sort()
list4.sort()

plt.figure(figsize=(16, 3))
plt.subplot(1,4,1)
plt.hist(list1, color='blue')
plt.xlabel('Patyway size')
plt.ylabel('Number of pathway')
plt.title('Disease pathways\n(reactome)')

plt.subplot(1,4,2)
plt.hist(list2, color='blue')
plt.xlabel('Patyway size')
plt.ylabel('Number of pathway')
plt.title('Disease pathways\n(non-reactome)')

plt.subplot(1,4,3)
plt.hist(list3, color='blue')
plt.xlabel('Patyway size')
plt.ylabel('Number of pathway')
plt.title('Comorbidity pathways\n(reactome)')

plt.subplot(1,4,4)
plt.hist(list4, color='blue')
plt.xlabel('Patyway size')
plt.ylabel('Number of pathway')
plt.title('Comorbidity pathways\n(non-reactome)')

# plt.show()
plt.savefig('a1.pdf', bbox_inches='tight')