import os


for flag in ['with', 'without']:
    if flag == 'with':
        path1 = '../genome/ldsc_file/ldsc_rg_with_HLA/'
        path2 = '../genome/ldsc_file/diseasePair_rg.txt'
    elif flag == 'without':
        path1 = '../genome/ldsc_file/ldsc_rg_without_HLA/'
        path2 = '../genome/ldsc_file/diseasePair_rg_rmhla.txt'

    file_list = os.listdir(path1)
    genetic_correlation = dict()
    for each in file_list:
        if not each.endswith('.log'):
            continue
        str1, str2 = each.replace('.log', '').split('~')
        key1 = str1.split('_')[-1]
        key2 = str2.split('_')[-1]
        if (len(key1) != 3) | (len(key2) != 3):
            continue
        with open(path1 + each, 'r') as infile:
            for line in infile:
                if 'Genetic Correlation:' in line:
                    str2 = line.split(':')
                    str3 = str2[1]
                    str4 = str3.split('(')[0]
                    rg = str4.replace(' ', '')
                    if rg == 'nan':
                        continue
                if 'P: ' in line:
                    p_value = line.strip('\r\n').split(' ')[1]
                    if p_value == 'nan':
                        continue
                    if key1 < key2:
                        genetic_correlation[key1 + '~' + key2] = rg + '\t' + p_value
                    else:
                        genetic_correlation[key2 + '~' + key1] = rg + '\t' + p_value
            infile.close()

    with open(path2, 'w+') as outfile:
        outfile.write('diseasePair\trg\tp\n')
        for each in genetic_correlation:
            outfile.write(each + '\t' + genetic_correlation[each] + '\n')
        outfile.close()