import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import simulation_test


if __name__ == '__main__':
    alleles_count_ = [10, 50, 100, 500, 1000]
    # alleles_count_ = [10, 50, 100, 500, 1000]
    population_amount_ = 10000000
    interval_for_alpha = 0.125

    alpha_values = np.arange(start=0.0, stop=1.0 + interval_for_alpha, step=interval_for_alpha)
    plot_means = []
    plot_stds = []

    index = [f'{alpha}' for alpha in alpha_values]
    columns = alleles_count_
    df_initial_times = pd.DataFrame(columns=columns, index=index, dtype=float)
    df_iterations_times = pd.DataFrame(columns=columns, index=index, dtype=float)

    for i, alpha in enumerate(alpha_values):
        for j, alleles in enumerate(alleles_count_):
            elapsed_initial_time, elapsed_iterations_time = simulation_test.perform_experiment(alleles, population_amount_, alpha, experiment_num=0)
            df_initial_times.iloc[i, j] = elapsed_initial_time
            df_iterations_times.iloc[i, j] = elapsed_iterations_time

    df_initial_times.to_csv(f'data/running_times/alleles_alphas_initial_times_population={population_amount_}')
    df_iterations_times.to_csv(f'data/running_times/alleles_alphas_iterations_times_population={population_amount_}')

    sns.heatmap(df_initial_times, annot=True, cmap=plt.cm.CMRmap_r)
    plt.xlabel('Alleles')
    plt.ylabel('Alpha values')
    plt.title(f'Initial running times')
 #   plt.savefig(f'plots/running_times/alleles_alphas_initial_times_population={population_amount_}',
 #               pad_inches=0.2, bbox_inches="tight")
    plt.show()

    sns.heatmap(df_iterations_times, annot=True, cmap=plt.cm.CMRmap_r)
    plt.xlabel('Alleles')
    plt.ylabel('Alpha values')
    plt.title(f'Iterations running times')
 #   plt.savefig(f'plots/running_times/alleles_alphas_iterations_times_population={population_amount_}',
 #               pad_inches=0.2, bbox_inches="tight")
    plt.show()