
all_variant = set()
with open('../genome/varaint_info.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        all_variant.add(line.split('\t')[0])
    infile.close()


snp_ld = dict()
with open('../LD_Calculation/1000_genome_download/1000genomes_from_plink2_grch37/1000genomes.ld', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = line.split()
        snp1 = str2[2]
        snp2 = str2[5]
        r2 = str2[6]
        if float(r2) <= 0.8:
            continue
        if snp1 in snp_ld:
            snp_ld[snp1].add(snp2)
        else:
            snp_ld[snp1] = set()
            snp_ld[snp1].add(snp2)
        if snp2 in snp_ld:
            snp_ld[snp2].add(snp1)
        else:
            snp_ld[snp2] = set()
            snp_ld[snp2].add(snp1)
    infile.close()


# disease_sig_snp
disease_snp_sig = dict()
with open('../genome/disease_snp_sig.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease = str2[0]
        snp = str2[1]
        if snp == '':
            continue
        disease_snp_sig[disease] = set(snp.split(';'))
    infile.close()

# disease_snp_sig_ld
disease_snp_sig_ld = dict()
for each in disease_snp_sig:
    disease = each
    snp_set = disease_snp_sig[disease]
    if '' in snp_set:
        snp_set.remove('')
    print(len(snp_set))
    set1 = snp_set
    for each1 in snp_set:
        if each1 in snp_ld:
            set1 = set1 | snp_ld[each1]
    print(len(set1))
    disease_snp_sig_ld[disease] = set1

# output
with open('../genome/disease_snp_sig_ld.txt', 'w+') as outfile:
    outfile.write('disease' + '\t' + 'snp' + '\n')
    for each in disease_snp_sig_ld:
        set1 = disease_snp_sig_ld[each] & all_variant
        outfile.write(each + '\t' + ';'.join(set1) + '\n')
    outfile.close()


# ----------------------------------- without hla --------------------------------- #

all_variant = set()
with open('../genome/varaint_info_rmhla.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        all_variant.add(line.split('\t')[0])
    infile.close()

# disease_sig_snp_rmhla
disease_snp_sig_rmhla = dict()
with open('../genome/disease_snp_sig_rmhla.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease = str2[0]
        snp = str2[1]
        if snp == '':
            continue
        disease_snp_sig_rmhla[disease] = set(snp.split(';'))
    infile.close()

# disease_snp_sig_rmhla_ld
disease_snp_sig_ld_rmhla = dict()
for each in disease_snp_sig_rmhla:
    disease = each
    snp_set = disease_snp_sig_rmhla[disease]
    if '' in snp_set:
        snp_set.remove('')
    print(len(snp_set))
    set1 = snp_set
    for each1 in snp_set:
        if each1 in snp_ld:
            set1 = set1 | snp_ld[each1]
    print(len(set1))
    disease_snp_sig_ld_rmhla[disease] = set1

# output
with open('../genome/disease_snp_sig_ld_rmhla.txt', 'w+') as outfile:
    outfile.write('disease' + '\t' + 'snp' + '\n')
    for each in disease_snp_sig_ld_rmhla:
        set1 = disease_snp_sig_ld_rmhla[each] & all_variant
        outfile.write(each + '\t' + ';'.join(set1) + '\n')
    outfile.close()