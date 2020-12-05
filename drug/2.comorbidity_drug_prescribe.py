import pandas as pd


gene_id_symbol = dict()
df = pd.read_table('../genome/grch37_download/genome_assemblies/gene_assembly_grch37.txt')
list1 = df.values.tolist()
for each in list1:
    id = str(each[0])
    symbol = str.upper(each[1])
    gene_id_symbol[id] = symbol


df = pd.read_table('../overlap/drug/ukb_drug_2_drugbank.txt')
list1 = df.values.tolist()
drug_durgbank_ukb_map = dict()
for each in list1:
    ukb_name = each[1]
    drugbank_id = each[2]
    if drugbank_id not in drug_durgbank_ukb_map:
        drug_durgbank_ukb_map[drugbank_id] = set()
    drug_durgbank_ukb_map[drugbank_id].add(ukb_name)


comorbidity_eat_drug_frequency = dict()
df = pd.read_csv('../overlap/drug/comorbidity_drug_score.csv')
list1 = df[['disease1', 'disease2', 'drug', 'frequency']].values.tolist()
for each in list1:
    code1 = each[0]
    code2 = each[1]
    drug = each[2]
    frequency = each[3]
    comorbidity_eat_drug_frequency[(code1, code2, drug)] = frequency


result = list()
comorbidity_drug_overlap = dict()
df = pd.read_csv('a.csv')
list1 = df.values.tolist()
for each in list1:
    code1, code2, gene, drug, drug_name, indication, treat1, treat2 = each
    symbol = gene_id_symbol[str(gene)]
    if drug in drug_durgbank_ukb_map:
        ukb_name = drug_durgbank_ukb_map[drug]
    else:
        result.append([code1, code2, symbol, drug, drug_name, indication, treat1, treat2, '', ''])
        continue
    for each1 in ukb_name:
        if (code1, code2, each1) in comorbidity_eat_drug_frequency:
            frequency = comorbidity_eat_drug_frequency[(code1, code2, each1)]
            result.append([code1, code2, symbol, drug, drug_name, indication, treat1, treat2, each1, frequency])
        else:
            result.append([code1, code2, symbol, drug, drug_name, indication, treat1, treat2, each1, ''])


df = pd.DataFrame(result, columns=['code1', 'code2', 'gene', 'drugbank_id', 'drugbank_name', 'indication', 'treat disease1', 'treat disease2', 'ukb_name', 'frequency'])

df.to_csv('b.csv', index=False)