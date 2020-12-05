import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np



df = pd.read_csv('../phenotype/comorbidity_filter.csv')
list1 = df[['disease1', 'disease2', 'RR']].values.tolist()

comorbidity_rr = dict()
for each in list1:
    code1 = each[0]
    code2 = each[1]
    rr = each[2]
    comorbidity_rr[(code1, code2)] = rr

comorbidity_rg = dict()
df = pd.read_table('../overlap/comorbidity_rg.txt')
list1 = df[['code1', 'code2', 'rg']].values.tolist()
for each in list1:
    code1 = each[0]
    code2 = each[1]
    rg = each[2]
    comorbidity_rg[(code1, code2)] = rg


rr_list = list()
rg_list = list()
for each in comorbidity_rr:
    rr = comorbidity_rr[each]
    if rr > 30:
        continue
    if each not in comorbidity_rg:
        continue
    rg = comorbidity_rg[each]
    rr_list.append(rr)
    rg_list.append(rg)

df = pd.DataFrame([rg_list, rr_list])
df = df.T
df.columns = ['rg', 'RR']
df1 = df[df['RR']>30]
print(df1.shape[0])
df2 = df.sort_values(by='rg')


[r, p] = stats.pearsonr(df2['rg'], df2['RR'])
print([r, p])


ax = sns.regplot(x="rg", y='RR', data=df2, x_jitter=.1, label=None,
                 marker='o', scatter_kws={'s': 3, 'alpha': 0.8}, color='dimgray', x_estimator=np.mean)
plt.text(0.1, 28, 'pearsonr=0.39', fontsize=16)
plt.text(0.1, 25, 'p_value=1.9e-72', fontsize=16)
plt.xlabel('rg', fontsize=16)
plt.ylabel('RR', fontsize=16)
ax.tick_params(labelsize=16)
plt.savefig('d.pdf', bbox_inches='tight')
# plt.show()