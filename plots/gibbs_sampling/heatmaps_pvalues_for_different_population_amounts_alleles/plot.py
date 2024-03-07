import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    # df_p_values_alpha_990 = pd.read_csv('df_p_values_alpha_990.csv', index_col=0)
    # df_p_values_alpha_1000 = pd.read_csv('df_p_values_alpha_1000.csv', index_col=0)
    df_p_values_small_population_alpha_920 = pd.read_csv('df_p_values_small_population_alpha_920.csv', index_col=0)
    df_p_values_small_population_alpha_1000 = pd.read_csv('df_p_values_small_population_alpha_1000.csv', index_col=0)

    # plot p_values
    # plt.figure(figsize=(10, 8))
    x_label = 'Population amounts'
    y_label = 'Alleles amounts'

    # # large population
    # plt.figure(figsize=(9, 6))
    # ax = sns.heatmap(df_p_values_alpha_990, annot=True, cmap='gray_r', annot_kws={"size": 14})
    # cax = ax.figure.axes[-1]
    # cax.tick_params(labelsize=14)
    # _, xlabels = plt.xticks()
    # ax.set_xticklabels(xlabels, size=14)
    # _, ylabels = plt.yticks()
    # ax.set_yticklabels(ylabels, size=14)
    #
    # plt.xlabel(x_label, fontsize=18)
    # plt.ylabel(y_label, fontsize=18)
    # # plt.title('P-values for alpha = 0.93')
    # plt.savefig('df_p_values_alpha_990', pad_inches=0.2, bbox_inches="tight")
    # plt.show()
    #
    # # plot times
    # plt.figure(figsize=(9, 6))
    # ax = sns.heatmap(df_p_values_alpha_1000, annot=True, cmap='gray_r', annot_kws={"size": 14})
    # cax = ax.figure.axes[-1]
    # cax.tick_params(labelsize=14)
    # _, xlabels = plt.xticks()
    # ax.set_xticklabels(xlabels, size=14)
    # _, ylabels = plt.yticks()
    # ax.set_yticklabels(ylabels, size=14)
    #
    # plt.xlabel(x_label, fontsize=18)
    # plt.ylabel(y_label, fontsize=18)
    # # plt.title('P-values for alpha = 1.0')
    # plt.savefig('df_p_values_alpha_1000', pad_inches=0.2, bbox_inches="tight")
    # plt.show()

    # small population
    plt.figure(figsize=(9, 6))
    ax = sns.heatmap(df_p_values_small_population_alpha_920, annot=True, cmap='gray_r', annot_kws={"size": 14})
    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=14)
    _, xlabels = plt.xticks()
    ax.set_xticklabels(xlabels, size=14)
    _, ylabels = plt.yticks()
    ax.set_yticklabels(ylabels, size=14)

    plt.xlabel(x_label, fontsize=18)
    plt.ylabel(y_label, fontsize=18)
    # plt.title('P-values for alpha = 0.93')
    plt.savefig('df_p_values_small_population_alpha_920', pad_inches=0.2, bbox_inches="tight")
    plt.show()

    # plot times
    plt.figure(figsize=(9, 6))
    ax = sns.heatmap(df_p_values_small_population_alpha_1000, annot=True, cmap='gray_r', annot_kws={"size": 14})
    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=14)
    _, xlabels = plt.xticks()
    ax.set_xticklabels(xlabels, size=14)
    _, ylabels = plt.yticks()
    ax.set_yticklabels(ylabels, size=14)

    plt.xlabel(x_label, fontsize=18)
    plt.ylabel(y_label, fontsize=18)
    # plt.title('P-values for alpha = 1.0')
    plt.savefig('df_p_values_small_population_alpha_1000', pad_inches=0.2, bbox_inches="tight")
    plt.show()