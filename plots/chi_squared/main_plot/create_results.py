import numpy as np
import pandas as pd
import csv
import os
import warnings
from models_simulated import chi_squared_with_csv
from models_simulated import chi_squared_second_attempt
from hwetests import asta


def save_results(alleles_amount, population_size, num_of_experiments, alpha_vals, uncertainty_vals):
    df_means_old = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
    df_means_correction = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
    df_means_asta = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)

    df_stds_old = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
    df_stds_correction = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
    df_stds_asta = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)

    # file_path = 'current_data.csv'

    for i in range(len(alpha_vals)):
        for j in range(len(uncertainty_vals)):
            results_old_list = []
            results_correction_list = []
            results_asta_list = []

            alleles_probabilities = chi_squared_second_attempt.generate_alleles_probabilities(
                alleles_count=alleles_amount)

            for experiment in range(num_of_experiments):
                # if experiment % 100 == 0:
                print(
                    f'experiment: {experiment}, population size: {population_size}, alleles amount: {alleles_amount}, alpha: {alpha_vals[i]}, '
                    f'uncertainty: {uncertainty_vals[j]}')
                print('generating data')

                data = chi_squared_second_attempt.generate_data(alleles_count=alleles_amount,
                                                                population_amount=population_size,
                                                                alleles_probabilities=alleles_probabilities,
                                                                alpha_val=alpha_vals[i],
                                                                uncertainty_val=uncertainty_vals[j])
                # print(len(data))
                print('running test')
                result_old, result_new, dof_old, dof_new = chi_squared_second_attempt.run_experiment(data=data,
                                                                                                     cutoff_value=2.0)
                print(f'alleles amount: {alleles_amount}, population: {population_size}')
                print(f'p_value_old: {result_old}, dof_old: {dof_old}')
                print('--------------')
                print(f'p_value_corrected: {result_new}, dof_corrected: {dof_new}')
                # chi_squared_second_attempt.plot_variance_vs_corrected_variance(data=data)
                with open('data.csv', 'w', newline='') as file:
                    writer = csv.writer(file)
                    for line in data:
                        writer.writerow(line)
                asta.full_algorithm(file_path='data.csv',
                                    cutoff_value=2.0,
                                    should_save_csv='simulation_data')
                # p_val, _, _ = asta.full_algorithm(file_path='data.csv',
                #                                   cutoff_value=2.0)
                result_asta = 0
                # print(f'')
                results_old_list.append(result_old)
                results_correction_list.append(result_new)
                results_asta_list.append(result_asta)
            # get mean and std
            mean_old = np.mean(results_old_list)
            std_old = np.std(results_old_list) / np.sqrt(len(results_old_list))

            mean_correction = np.mean(results_correction_list)
            std_correction = np.std(results_correction_list) / np.sqrt(len(results_correction_list))

            mean_asta = np.mean(results_asta_list)
            std_asta = np.std(results_asta_list) / np.sqrt(len(results_asta_list))

            df_means_old.iloc[i, j] = mean_old
            df_stds_old.iloc[i, j] = std_old

            df_means_correction.iloc[i, j] = mean_correction
            df_stds_correction.iloc[i, j] = std_correction

            df_means_asta.iloc[i, j] = mean_asta
            df_stds_asta.iloc[i, j] = std_asta

    # path to this directory
    current_path = os.getcwd()

    # -> make directory: data -> real_data -> levels
    path_to_save = os.path.join(current_path, f'population_{population_size}_alleles_{alleles_amount}_test')
    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save)

    # save dataframes
    df_means_old.to_csv(f'{path_to_save}/means_old.csv')
    df_stds_old.to_csv(f'{path_to_save}/stds_old.csv')

    df_means_correction.to_csv(f'{path_to_save}/means_correction.csv')
    df_stds_correction.to_csv(f'{path_to_save}/stds_correction.csv')

    df_means_asta.to_csv(f'{path_to_save}/means_asta.csv')
    df_stds_asta.to_csv(f'{path_to_save}/stds_asta.csv')


if __name__ == '__main__':
    warnings.filterwarnings("error")

    experiments_amount = 1
    # alleles_amounts = [50, 100, 200, 500]  # 2
    # population_sizes = [50000, 100000, 100000, 100000]  # 35
    alleles_amounts = [500]  # 2
    population_sizes = [100000]  # 35

    interval_for_alpha = 0.04  # 0.02
    # uncertainty = 0.2
    uncertainty_values = [0.0, 0.1, 0.2, 0.4]

    alpha_values = np.arange(start=0.92, stop=1 + interval_for_alpha,
                             step=interval_for_alpha)  # start, stop, step
    alpha_values = np.array([round(alpha, 2) for alpha in alpha_values])

    alpha_values = [1.0]
    uncertainty_values = [0.4]

    for i_ in range(len(alleles_amounts)):
        alleles_amount_ = alleles_amounts[i_]
        population_size_ = population_sizes[i_]
        save_results(alleles_amount=alleles_amount_,
                     population_size=population_size_,
                     num_of_experiments=experiments_amount,
                     alpha_vals=alpha_values,
                     uncertainty_vals=uncertainty_values)
