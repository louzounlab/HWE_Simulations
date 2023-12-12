import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def change_index(index):
    vocab = {}
    for i in index:
        prefix = i / 1000000
        if prefix.is_integer():
            prefix = int(prefix)
        vocab[i] = str(prefix) + 'M'
        if str(prefix) == '1000':
            vocab[i] = '1B'
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
    plt.ylabel('Population amounts', fontsize=18)
    # plt.title('P-values')
    # plt.title('C', weight='bold', fontsize=20)
    plt.savefig('heatmap_populations_alphas.pdf', format='pdf', bbox_inches="tight")
    plt.show()