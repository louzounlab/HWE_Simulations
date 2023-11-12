import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    sns.set_style('white')
    # sns.set(font_scale=1.3)
    plt.rcParams["font.family"] = "Arial"

    df_times = pd.read_csv('df_times.csv', index_col=0)

    alleles = df_times.index.to_list()
    times_1m = df_times['1000000'].to_list()
    times_10m = df_times['10000000'].to_list()
    times_100m = df_times['100000000'].to_list()

    plt.figure(figsize=(9, 6))
    plt.plot(alleles, times_1m, '-o', label='1M population size', color='slategrey')
    plt.plot(alleles, times_10m, '-o', label='10M population size', color='limegreen')
    plt.plot(alleles, times_100m, '-o', label='100M population size', color='deeppink')

    plt.xlabel('Alleles amounts', fontsize=18)
    plt.ylabel('Elapsed time in seconds', fontsize=18)
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    # min_ = min(chi_squared_results, gibbs_sampling_results)
    # max_ = max(chi_squared_results, gibbs_sampling_results)
    # print(gibbs_sampling_results)
    # print(chi_squared_results)
    # plt.xlim(0.0, 1.0)
    # plt.ylim(0.0, 1.0)
    # plt.title('D', weight='bold', fontsize=20)
    plt.legend(fontsize=14)
    plt.savefig('df_times.png', pad_inches=0.2, bbox_inches="tight")
    plt.show()