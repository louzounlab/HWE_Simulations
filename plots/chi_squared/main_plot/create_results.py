import numpy as np
import pandas as pd
import os
import warnings
from models_simulated import chi_squared_with_csv
from models_simulated import chi_squared_second_attempt


def save_results(alleles_amount, population_size, num_of_experiments, alpha_vals, uncertainty_vals):
    df_means_old = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
    df_means_correction = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)

    df_stds_old = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
    df_stds_correction = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)

    # file_path = 'current_data.csv'

    for i in range(len(alpha_vals)):
        for j in range(len(uncertainty_vals)):
            results_old_list = []
            results_correction_list = []

            alleles_probabilities = chi_squared_second_attempt.generate_alleles_probabilities(alleles_count=alleles_amount)

            for experiment in range(num_of_experiments):
                # if experiment % 100 == 0:
                print(f'experiment: {experiment}, population size: {population_size}, alleles amount: {alleles_amount}, alpha: {alpha_vals[i]}, '
                      f'uncertainty: {uncertainty_vals[j]}')
                print('generating data')
                data = chi_squared_second_attempt.generate_data(alleles_count=alleles_amount,
                                                                population_amount=population_size,
                                                                alleles_probabilities=alleles_probabilities,
                                                                alpha_val=alpha_vals[i],
                                                                uncertainty_val=uncertainty_vals[j])
                # print(len(data))
                print('running test')
                result_old, result_new = chi_squared_second_attempt.run_experiment(data=data,
                                                                                   cutoff_value=2.0)
                results_old_list.append(result_old)
                results_correction_list.append(result_new)
            # get mean and std
            mean_old = np.mean(results_old_list)
            std_old = np.std(results_old_list) / np.sqrt(len(results_old_list))

            mean_correction = np.mean(results_correction_list)
            std_correction = np.std(results_correction_list) / np.sqrt(len(results_correction_list))

            df_means_old.iloc[i, j] = mean_old
            df_stds_old.iloc[i, j] = std_old

            df_means_correction.iloc[i, j] = mean_correction
            df_stds_correction.iloc[i, j] = std_correction

    # path to this directory
    current_path = os.getcwd()

    # -> make directory: data -> real_data -> levels
    path_to_save = os.path.join(current_path, f'population_{population_size}_alleles_{alleles_amount}')
    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save)

    # save dataframes
    df_means_old.to_csv(f'{path_to_save}/means_old.csv')
    df_stds_old.to_csv(f'{path_to_save}/stds_old.csv')

    df_means_correction.to_csv(f'{path_to_save}/means_correction.csv')
    df_stds_correction.to_csv(f'{path_to_save}/stds_correction.csv')


if __name__ == '__main__':
    warnings.filterwarnings("error")

    experiments_amount = 1000
    alleles_amounts = [50, 100, 200, 500]  # 2
    population_sizes = [50000, 100000, 100000, 100000]  # 35

    interval_for_alpha = 0.04  # 0.02
    # uncertainty = 0.2
    uncertainty_values = [0.0, 0.1, 0.2, 0.4]

    alpha_values = np.arange(start=0.92, stop=1 + interval_for_alpha,
                             step=interval_for_alpha)  # start, stop, step
    alpha_values = np.array([round(alpha, 2) for alpha in alpha_values])

    for i in range(len(alleles_amounts)):
        alleles_amount_ = alleles_amounts[i]
        population_size_ = population_sizes[i]
        save_results(alleles_amount=alleles_amount_,
                     population_size=population_size_,
                     num_of_experiments=experiments_amount,
                     alpha_vals=alpha_values,
                     uncertainty_vals=uncertainty_values)
