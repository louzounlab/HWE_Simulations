import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import single_experiment_with_certainty
import time

if __name__ == '__main__':
    alleles_nums = np.array([10, 250])
    population_nums = np.array([1000000, 10000000, 100000000])
    experiments_amount = 10

    alpha_vals = [0.0, 0.5]
    for alpha_val in alpha_vals:
        # prepare dataframes
        index = [f'{alleles}' for alleles in alleles_nums]
        columns = [f'{population}' for population in population_nums]
        initial_times_df = pd.DataFrame(columns=columns, index=index, dtype=float)
        iterations_times_df = pd.DataFrame(columns=columns, index=index, dtype=float)
        positive_results_df = pd.DataFrame(columns=columns, index=index, dtype=float)

        for i, alleles in enumerate(alleles_nums):
            for j, population in enumerate(population_nums):
                # prepare the simulation data
                population_amount_calculated, alleles_probabilities, probabilities, observed, elapsed_initial_time = \
                    single_experiment_with_certainty.prepare_experiment_data(alleles_count=alleles,
                                                                             population_amount=population,
                                                                             alpha_val=alpha_val)
                # with the simulation data, perform the simulation 50 times.
                # count the time it took to run all 50 simulations
                start_time = time.time()

                # we have observed matrix, in each experiment we change it, but we want the same matrix for all
                # experiments, therefore we will save a copy (observed_copy)
                # and in each new experiment we will just copy observed to observed_copy using copyto()
                # and use observed_copy in each experiment.

                observed_copy = np.copy(observed)

                results = []
                for experiment_num in range(experiments_amount):
                    result = \
                        single_experiment_with_certainty.perform_experiment(alleles_count=alleles,
                                                                            population_amount=population,
                                                                            alpha_val=alpha_val,
                                                                            experiment_num=experiment_num,
                                                                            population_amount_calculated=population_amount_calculated,
                                                                            alleles_probabilities=alleles_probabilities,
                                                                            probabilities=probabilities,
                                                                            observed=observed_copy)
                    results.append(result)
                    # initialize observed_copy, copy from observed matrix
                    np.copyto(observed_copy, observed)
                # we have a list for the results, and the amount of time it took to run 50 results.
                # take average value for the results list and save in dataframes.
                mean_result = np.mean(results)
                end_time = time.time()
                elapsed_iterations_time = end_time - start_time

                initial_times_df.iloc[i, j] = elapsed_initial_time
                iterations_times_df.iloc[i, j] = elapsed_iterations_time
                positive_results_df.iloc[i, j] = mean_result
        # finished saving data in dataframes, now save them as csvs.
        # convert alpha -> alpha * 100 for naming files with alpha (can't use a dot)
        alpha_val_format = int(alpha_val * 100)
        initial_times_df.to_csv(f'data/simulation_data/gibbs_sampling_data/test_initial_times_alpha={alpha_val_format}')
        iterations_times_df.to_csv(f'data/simulation_data/gibbs_sampling_data/test_simulation_times_alpha={alpha_val_format}')
        positive_results_df.to_csv(f'data/simulation_data/gibbs_sampling_data/test_results_alpha={alpha_val_format}')