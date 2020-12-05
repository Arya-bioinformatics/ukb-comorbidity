import gzip
import os


root_dir = '../geneAtlas/all_variant/'
file_list = os.listdir(root_dir)
HLA_region = 'chr6:29691116–33054976'

list1 = list()
list2 = list()
set1 = set()
list1.append('snp' + '\t' + 'chr' + '\t' + 'position' + '\t' + 'A1' + '\t' + 'A2' + '\t' + 'MAF' + '\t' + 'HWE_P' + '\t' + 'imputed_score')
list2.append('snp' + '\t' + 'chr' + '\t' + 'position' + '\t' + 'A1' + '\t' + 'A2' + '\t' + 'MAF' + '\t' + 'HWE_P' + '\t' + 'imputed_score')
for each in file_list:
    if not each.endswith('csv.gz'):
        continue
    elif 'genotyped' in each:
        chr = each.split('.')[2].replace('chr', '')
        f = gzip.open(root_dir + each, 'rb')
        file_content = str(f.read(), encoding = "utf8")
        all_line = file_content.split('\n')
        for line in all_line[1:]:
            str1 = line.split(' ')
            if len(str1) == 1:
                continue
            rsid = str1[0]
            if not rsid.startswith('rs'):
                continue
            if rsid in set1:
                continue
            set1.add(rsid)
            position = str1[1]
            a1 = str1[2]
            a2 = str1[3]
            maf = str1[4]
            HWE = str1[5]
            if float(maf) < 0.01:
                continue
            if float(HWE) < float('1e-50'):
                continue
            list1.append(rsid + '\t' + chr + '\t' + position + '\t' + a1 + '\t' + a2 + '\t' + maf + '\t' + HWE + '\t' + '')
            if (chr=='6') & (int(position)>=29691116) & (int(position)<=33054976):
                continue
            else:
                list2.append(rsid + '\t' + chr + '\t' + position + '\t' + a1 + '\t' + a2 + '\t' + maf + '\t' + HWE + '\t' + '')
        f.close()
    elif 'imputed' in each:
        print(each)
        chr = each.split('.')[2].replace('chr', '')
        f = gzip.open(root_dir + each, 'rb')
        file_content = str(f.read(), encoding="utf8")
        all_line = file_content.split('\n')
        for line in all_line[1:]:
            str1 = line.split(' ')
            if len(str1) == 1:
                continue
            rsid = str1[0]
            if not rsid.startswith('rs'):
                continue
            if rsid in set1:
                continue
            set1.add(rsid)
            position = str1[1]
            a1 = str1[2]
            a2 = str1[3]
            maf = str1[4]
            HWE = str1[5]
            imputed_score = str1[6]
            if float(maf) < 0.01:
                continue
            if float(HWE) < float('1e-50'):
                continue
            if float(imputed_score) < 0.9:
                continue
            list1.append(rsid + '\t' + chr + '\t' + position + '\t' + a1 + '\t' + a2 + '\t' + maf + '\t' + HWE + '\t' + imputed_score)
            if (chr=='6') & (int(position)>=29691116) & (int(position)<=33054976):
                continue
            else:
                list2.append(rsid + '\t' + chr + '\t' + position + '\t' + a1 + '\t' + a2 + '\t' + maf + '\t' + HWE + '\t' + imputed_score)
        f.close()

with open('../genome/varaint_info.txt', 'w+') as outfile:
    outfile.write('\n'.join(list1))
    outfile.close()

with open('../genome/varaint_info_rmhla.txt', 'w+') as outfile:
    outfile.write('\n'.join(list2))
    outfile.close()