import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

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
    with open('../MimMiner_download/MimMiner_Exp_AC_T_TXCS_basedonACMESH_filt_RW.txt',
              'r') as infile:
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

    df = pd.read_csv('../EMR/phenotype/comorbidity_filter.csv')
    list1 = df[['disease1', 'disease2', 'RR']].values.tolist()

    comorbidity_rr = dict()
    for each in list1:
        code1 = each[0]
        code2 = each[1]
        rr = each[2]
        # if rr > 60:
        #     continue
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
        comorbidity_similarity.append([code1, code2, float(rr), float(s)])

    df = pd.DataFrame(comorbidity_similarity, columns=['disease1', 'disease2', 'rr', 'PheSim'])
    df1 = df[df['rr']>60]
    print(df1.shape[0])
    df.to_csv('a.csv', index=False)


    df1 = df[['rr', 'PheSim']]
    df1 = df1.sort_values(by='PheSim')


    [r, p] = stats.pearsonr(df1['rr'], df1['PheSim'])
    print([r, p])


    ax = sns.regplot(x="PheSim", y='rr', data=df1, x_jitter=.1, label=None,
                     marker='o', scatter_kws={'s': 5, 'alpha': 0.8}, color='dimgray', x_estimator=np.mean)
    plt.text(0.0, 60, 'pearsonr=0.27', fontsize=16)
    plt.text(0.0, 55, 'p_value=1.5e-8', fontsize=16)
    plt.xlabel('PheSim', fontsize=16)
    plt.ylabel('RR', fontsize=16)
    ax.tick_params(labelsize=16)
    plt.savefig('d.pdf', bbox_inches='tight')
    # plt.show()