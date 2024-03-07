import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def change_index(index):
    vocab = {}
    for i in index:
        prefix = i / 1000
        if prefix.is_integer():
            prefix = int(prefix)
            vocab[i] = str(prefix) + 'K'
        else:
            vocab[i] = str(i)
    return vocab


if __name__ == '__main__':
    sns.set_style('white')
    # sns.set(font_scale=1.3)
    plt.rcParams["font.family"] = "Arial"

    df_p_values = pd.read_csv('df_p_values.csv', index_col=0)
    df_p_values.rename(index=change_index(df_p_values.index), inplace=True)
    # df_times = pd.read_csv('df_times.csv', index_col=0)

    # plot p_values
    plt.figure(figsize=(9, 6))
    ax = sns.heatmap(df_p_values, annot=True, fmt='.3g', cmap='gray_r', annot_kws={"size": 14})
    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=14)

    _, xlabels = plt.xticks()
    ax.set_xticklabels(xlabels, size=14)

    _, ylabels = plt.yticks()
    ax.set_yticklabels(ylabels, size=14)

    plt.xlabel('Alpha values', fontsize=18)
    plt.ylabel('Alleles amounts', fontsize=18)
    # plt.title('P-values')
    # plt.title('C', weight='bold', fontsize=20)
    plt.savefig('df_p_values.png', pad_inches=0.2, bbox_inches="tight")
    plt.show()

    # plot times
    # sns.heatmap(df_times, annot=True, cmap='Blues')
    # plt.xlabel('Alpha values')
    # plt.ylabel('Alleles amounts')
    # plt.title('Elapsed time in seconds')
    # plt.savefig('df_times', pad_inches=0.2, bbox_inches="tight")
    # plt.show()
    # alleles = df_times.index.to_list()
    # times = df_times['1.0'].to_list()
    # plt.figure(figsize=(9, 6))
    # plt.plot(alleles, times, '-o')
    # # plt.scatter(gibbs_sampling_results, gibbs_sampling_results, s=4, color='deeppink', label='Gibbs Sampling')  # hotpink
    # # plt.title('')
    # # plt.xlabel('Snipping indices (normalized between 0 and 1)')
    # # plt.ylabel('p-value')
    # plt.xlabel('Alleles amounts', fontsize=16)
    # plt.ylabel('Elapsed time in seconds', fontsize=16)
    # ax = plt.gca()
    # ax.tick_params(axis='x', labelsize=12)
    # ax.tick_params(axis='y', labelsize=12)
    # # min_ = min(chi_squared_results, gibbs_sampling_results)
    # # max_ = max(chi_squared_results, gibbs_sampling_results)
    # # print(gibbs_sampling_results)
    # # print(chi_squared_results)
    # # plt.xlim(0.0, 1.0)
    # # plt.ylim(0.0, 1.0)
    # # plt.title('D', weight='bold', fontsize=20)
    # plt.savefig('df_times.png', pad_inches=0.2, bbox_inches="tight")
    # plt.show()