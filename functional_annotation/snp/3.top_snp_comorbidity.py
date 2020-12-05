import pandas as pd


for flag in ['with', 'without']:
    print('\n' + flag + '\n')
    if flag == 'with':
        path1 = '../overlap/comorbidity_snp.txt'
    if flag == 'without':
        path1 = '../overlap/comorbidity_snp_rmhla.txt'

    comorbidity_snp = dict()
    snp_comorbidity_count = dict()
    disease_code_name = dict()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            set1 = set(str2[4].split(';'))
            disease_code_name[str2[0]] = str2[2]
            disease_code_name[str2[1]] = str2[3]
            comorbidity_snp[(str2[0], str2[1])] = set1
            for each in set1:
                if each not in snp_comorbidity_count:
                    snp_comorbidity_count[each] = 0
                snp_comorbidity_count[each] += 1
        infile.close()

    list1 = list()
    for each in snp_comorbidity_count:
        list1.append([each, snp_comorbidity_count[each]])
    df = pd.DataFrame(list1, columns=['snp', 'count'])
    df = df.sort_values(by='count', ascending=False)

    top_snp = set()
    list1 = df.values.tolist()
    for each in list1:
        top_snp.add(each[0])
        if len(top_snp) == 10:
            break

    involved_comorbidity = set()
    for each in comorbidity_snp:
        name1 = disease_code_name[each[0]]
        name2 = disease_code_name[each[1]]
        if len(top_snp & comorbidity_snp[each]) != 0:
            involved_comorbidity.add(name1 + ' ------- ' + name2)

    for each in involved_comorbidity:
        print(each)
    print('total snp interpreted comorbidity: ' + str(len(comorbidity_snp)))
    print('total top snp interpreted comorbidity: ' + str(len(involved_comorbidity)))
    print('ratio: ' + str(len(involved_comorbidity)/len(comorbidity_snp)))