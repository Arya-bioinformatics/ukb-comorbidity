

for flag in ['with', 'without']:
    if flag == 'with':
        path1 = '../genome/varaint_info.txt'
        path2 = '../functional_annotation/annovar_input.avinput'
        path3 = '../functional_annotation/polyphen_input.txt'
        path4 = '..functional_annotation/sift_input.txt'
    if flag == 'without':
        path1 = '../genome/varaint_info_rmhla.txt'
        path2 = '../functional_annotation/annovar_input_rmhla.avinput'
        path3 = '../functional_annotation/polyphen_input_rmhla.txt'
        path4 = '../functional_annotation/sift_input.txt'

    variant_info = dict()
    set1 = set()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            rsid = str2[0]
            chromosome = str2[1]
            position = str2[2]
            a1 = str2[3]
            a2 = str2[4]
            variant_info[rsid] = [chromosome, position, position, a1, a2]
            set1.add(rsid)
        infile.close()

    # annovar input
    with open(path2, 'w+') as outfile:
        for each in variant_info:
            outfile.write('\t'.join(variant_info[each]) + '\t' + each + '\n')
        outfile.close()

    # polyphen 2 input
    with open(path3, 'w+') as outfile:
        outfile.write('\n'.join(set1))
        outfile.close()

    # sift input
    with open(path4, 'w+') as outfile:
        outfile.write('\n'.join(set1))
        outfile.close()