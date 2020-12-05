import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.stats import pearsonr
import seaborn as sns
import numpy as np

if __name__ == '__main__':

    disease_prevalence = dict()
    df = pd.read_excel('../Supplementary Data 1.xlsx')
    list1 = df.values.tolist()
    for each in list1:
        disease_prevalence[each[0]] = each[3]/410309

    df = pd.read_csv('../phenotype/comorbidity_filter.csv')
    list1 = df.values.tolist()
    disease_comor_count = dict()
    for each in list1:
        code1 = each[0]
        code2 = each[1]
        if code1 not in disease_comor_count:
            disease_comor_count[code1] = 0
        if code2 not in disease_comor_count:
            disease_comor_count[code2] = 0
        disease_comor_count[code1] += 1
        disease_comor_count[code2] += 1

    df = pd.DataFrame.from_dict(disease_comor_count, orient='index', columns=['count'])
    df = df.sort_values(by='count', ascending=False)
    df.to_csv('a.csv')

    list1 = list()
    for each in disease_comor_count:
        count = disease_comor_count[each]
        list1.append(math.log10(count))

    i = 0
    for each in list1:
        if each > 100:
            i += 1

    print('disease with more than 100 comorbidities: ' + str(i) + '\n')

    df = pd.DataFrame(list1, columns=['count'])
    df.to_csv('a.csv', index=False)

    plt.figure(figsize=(6, 4))
    plt.hist(list1, bins=10, normed=False, alpha=0.5, histtype='bar',
             color='steelblue', edgecolor='black', linewidth=0.5)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Number of comorbidities (log10)', fontsize=15)
    plt.ylabel('Number of diseases', fontsize=15)
    plt.savefig('a.png', bbox_inches='tight')

    list1 = list()
    for each in disease_comor_count:
        prevalence = disease_prevalence[each]
        count = disease_comor_count[each]
        list1.append([each, prevalence, count])
    df = pd.DataFrame(list1, columns=['disease', 'prevalence', 'count'])
    df = df.sort_values(by='prevalence')
    [coef, p] = pearsonr(df['prevalence'], df['count'])
    print([coef, p])

    plt.figure(figsize=(6, 4))
    ax = sns.regplot(x="prevalence", y='count', data=df, x_jitter=.1, label=None,
                     marker='o', scatter_kws={'s': 3, 'alpha': 0.8}, color='dimgray', x_estimator=np.mean)
    plt.text(0, 600, 'pearsonr=0.69', fontsize=16)
    plt.text(0, 550, 'p_value=3.3e-64', fontsize=16)
    plt.xlabel('prevalence', fontsize=16)
    plt.ylabel('comorbidity count', fontsize=16)
    ax.tick_params(labelsize=16)
    # plt.show()
    plt.savefig('a1.pdf', bbox_inches='tight')