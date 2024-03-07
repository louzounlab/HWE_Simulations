import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

if __name__ == '__main__':
    sns.set_style('white')
    plt.rcParams["font.family"] = "Arial"

    plt.figure(figsize=(9, 6))

    # Read results Dataframe
    df = pd.read_csv('results.csv', index_col=0)
    snipping_values = df['snipping_index'].to_list()
    gibbs_sampling_results = df['gibbs_sampling_p_value'].to_list()
    chi_squared_results = df['chi_squared_p_value'].to_list()
    snipping_amount = len(snipping_values)

    # Make plot
    snipping_indices = [snipping / (snipping_amount - 1) for snipping in range(snipping_amount)]
    plt.scatter(gibbs_sampling_results, chi_squared_results, color='black')
    # plt.scatter(gibbs_sampling_results, gibbs_sampling_results, s=4, color='deeppink', label='Gibbs Sampling')  # hotpink
    # plt.title('P-values')
    # plt.title('A', weight='bold', fontsize=20)
    # plt.xlabel('Snipping indices (normalized between 0 and 1)')
    # plt.ylabel('p-value')
    plt.xlabel('Gibbs sampling results', fontsize=18)
    plt.ylabel('Chi-Squared results', fontsize=18)
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    min_ = min(chi_squared_results, gibbs_sampling_results)
    max_ = max(chi_squared_results, gibbs_sampling_results)
    print(gibbs_sampling_results)
    print(chi_squared_results)
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.savefig('snipping_data_plot.pdf', format='pdf', bbox_inches="tight")
    plt.show()