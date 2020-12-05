import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


icd10_omim_map = dict()
with open('../name_map/all_omim_icd10_map.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        omim = str2[0]
        icd10 = str2[2][:3]
        icd10_omim_map[icd10] = omim
    infile.close()

# 按顺序读入相似性矩阵中的omim id
omimID_list = list()
with open('../MimMiner_download/MimMiner_Exp_AC_T_TXCS_basedonACMESH_filt_RW.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        omimID_list.append(str2[0])
    infile.close()

# 将相似性矩阵以字典形式存储
similarity_dict = dict()
with open('../MimMiner_download/MimMiner_Exp_AC_T_TXCS_basedonACMESH_filt_RW.txt', 'r') as infile:
    i = 0
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        omimID = str2[0]
        j = i
        for each in str2[j + 1:]:
            key = omimID + '*' + omimID_list[j]
            similarity_dict[key] = float(str2[j + 1])
            j = j + 1
        i = i + 1
    infile.close()

comorbidity_similarity = list()
x = 0
df = pd.read_table('../overlap/comorbidity_rg.txt')
list1 = df[['code1', 'code2', 'rg']].values.tolist()
for each in list1:
    code1 = each[0]
    code2 = each[1]
    if code1 not in icd10_omim_map:
        continue
    if code2 not in icd10_omim_map:
        continue
    omim_1 = icd10_omim_map[code1]
    omim_2 = icd10_omim_map[code2]
    key1 = omim_1 + '*' + omim_2
    key2 = omim_2 + '*' + omim_1
    if key1 in similarity_dict:
        s = similarity_dict[key1]
    elif key2 in similarity_dict:
        s = similarity_dict[key2]
    else:
        continue
    rg = each[2]
    if float(s) > 0.6:
        x += 1
        continue
    comorbidity_similarity.append([code1, code2, float(rg), float(s)])

df = pd.DataFrame(comorbidity_similarity, columns=['disease1', 'disease2', 'rg', 'PheSim'])
df.to_csv('a.csv', index=False)


df1 = df[['rg', 'PheSim']]
df1 = df1.sort_values(by='rg')


[r, p] = stats.pearsonr(df1['rg'], df1['PheSim'])
print([r, p])


ax = sns.regplot(x="rg", y='PheSim', data=df1, x_jitter=.1, label=None,
                 marker='o', scatter_kws={'s': 5, 'alpha': 0.8}, color='dimgray', x_estimator=np.mean)
plt.text(0.15, 0.43, 'pearsonr=0.22', fontsize=16)
plt.text(0.15, 0.40, 'p_value=0.03', fontsize=16)
plt.xlabel('rg', fontsize=16)
plt.ylabel('PheSim', fontsize=16)
ax.tick_params(labelsize=16)
# plt.show()
plt.savefig('c.pdf', bbox_inches='tight')