import numpy as np
import random

import pandas as pd
import seaborn as sns
import array as arr

import utils_with_certainty
import matplotlib.pyplot as plt
import time


# populations_probs of size alleles_count x alleles_count x population_amount
def perform_experiment(alleles_count, population_amount, alpha_val, experiment_num):
    # alleles_individuals = np.zeros(
    #    (population_amount, 2))  # matrix where every row is the alleles for a person
    print(f'alleles: {alleles_count}, population: {population_amount}')
    start_time = time.time()

    # amount of times we got O(i,j)=0 and needed to calculate again

    # probabilities {p(i)}
    alleles_probabilities = utils_with_certainty.calculate_alleles_probabilities(alleles_count)

    # checking how many unbalanced alleles we have
    # unbalanced_alleles = np.zeros(alleles_count)
    # unbalanced_alleles_counters = []

    # probabilities: {p(i,j)}. only upper triangle
    probabilities = utils_with_certainty.calc_p_ij(alleles_count, alleles_probabilities)
    # print(probabilities)

    # observed: {O(i,j) = Poisson()}. only upper triangle. here alpha affects the non-diagonal elements
    observed = utils_with_certainty.calc_O_ij(population_amount, alleles_count, alpha_val, alleles_probabilities,
                                              probabilities)
    print('calculated observed')
    print(observed)

    # calculate new amount of population based on the observations
    population_amount_calculated = 0
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            population_amount_calculated += observed[i, j]

    # # update dynamic alleles probabilities:
    for i in range(alleles_count):
        sum_current = 0
        for j in range(alleles_count):
            if i == j:
                sum_current += utils_with_certainty.get_O_ij(observed, i, j)
            else:
                sum_current += 0.5 * utils_with_certainty.get_O_ij(observed, i, j)
        alleles_probabilities[i] = sum_current / population_amount_calculated

    # [ln(p_0), ln(p_1),...,]
    # actually [1, delta_1, delta_2,...,delta_k]
    list_probabilities = []

    # calculate cdf. [(p_11, 1), (p_11+p_12, 2), ..., (p_11+...+p_kk, k)]
    cdf_dict = {}
    utils_with_certainty.calculate_cdf_dict(alleles_count, observed, cdf_dict)

    # add first ln(probability) as 0. We only care about the deltas anyway.
    list_probabilities.append(0)

    # calc couples
    first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
    second_couple = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                                   k=2)

    # pick two alleles using the cdf
    # (generating a uniform[0,1] r.v. and use it to generate r.v. with the given probabilities
    # using the Inverse cdf Theorem with binary search on the cdf list)
    # t
    allele_1 = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                              k=1)[0]
    # m
    allele_2 = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
                                                             left_alleles=second_couple)

    couple_from_cdf = [allele_1, allele_2]

    # every row represents a couple i,j and a value (either 1 or -1)
    # in the first iteration here we have 4 couples, but in the next iterations we will have 2 couples.
    couples = np.zeros(shape=(4, 3), dtype=int)
    couples[0, :] = [first_couple[0], first_couple[1], -1]
    couples[1, :] = [second_couple[0], second_couple[1], +1]
    couples[2, :] = [first_couple[1], allele_1, +1]
    couples[3, :] = [second_couple[1], allele_2, -1]

    # calculate new delta and modify observed
    utils_with_certainty.update_current_delta_probability(list_probabilities, observed, alleles_probabilities, couples,
                                                          population_amount_calculated, probabilities)

    end_time = time.time()
    elapsed_initial_time = end_time - start_time
    start_time = time.time()

    # now keep iterating to fill the list of deltas
    # iterations = 10000 * 10
    iterations = 100000
    for k in range(iterations):

        # unbalanced_alleles_counters.append(np.sum(np.abs(unbalanced_alleles)))

        # if k % 10000 == 0:
        #     cdf_dict = utils_with_certainty.calculate_cdf_dict(alleles_count, observed)

        # sums.append(sum_loops[0])
        print(f' loop {k} / {iterations}. experiment num:{experiment_num}. alleles: {alleles_count}, '
              f'population: {population_amount}. Alpha: {alpha_val}')

        # updating the couples
        first_couple = [first_couple[0], couple_from_cdf[1]]
        second_couple = [second_couple[0], couple_from_cdf[0]]

        # pick two alleles using the cdf (first we pick a number between 0 and 1 and then get element with the closest
        # probability)
        # x
        x = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                           k=1)[0]
        # y
        y = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
                                                          left_alleles=second_couple)

        couple_from_cdf = [x, y]

        # print(f'first_couple: [{first_couple[0], first_couple[1]}]')
        # print(f'second_couple: [{second_couple[0], second_couple[1]}]')
        # print(f'couple_from_cdf: [{couple_from_cdf[0], couple_from_cdf[1]}]')

        couples = np.zeros(shape=(2, 3), dtype=int)
        couples[0, :] = [first_couple[1], x, +1]
        couples[1, :] = [second_couple[1], y, -1]

        # calculate new delta
        # list_probabilities.append(utils_with_certainty.calc_current_delta_probability(observed, probabilities, couples,
        #                                                                               population_amount_calculated))
        utils_with_certainty.update_current_delta_probability(list_probabilities, observed, alleles_probabilities,
                                                              couples,
                                                              population_amount_calculated, probabilities)
        # utils_with_certainty.debug_observed(observed)
        if k == iterations - 1:
            print(f'final observed:')
            print(observed)

    # now we have the list of probabilities. check if 95% of the elements (sum of deltas) are bigger than 1.
    sum_current = 0
    bigger_counter = 0
    start_from = 10000

    values = []

    for delta in list_probabilities[:start_from]:
        sum_current += delta

    for delta in list_probabilities[start_from:]:
        sum_current += delta
        values.append(sum_current)
        if sum_current > list_probabilities[0]:
            bigger_counter += 1
    # return bigger_counter / (len(list_probabilities) - 1)
    # if (bigger_counter / len(list_probabilities) - 1) > 0.95:
    #     return 1
    # return 0
    # print(f'len: {len(values)}. values: {values}')

    plt.scatter(list(range(start_from, len(values) + start_from)), values, marker='.')
    plt.xlabel('Time')
    plt.ylabel('Approximated ln(p_i)')
    plt.title(f'Alpha: {alpha_val}. Alleles: {alleles_count}. Population: {population_amount}')
    plt.show()

    # return bigger_counter / (len(list_probabilities) - start_from)
    end_time = time.time()
    elapsed_iterations_time = end_time - start_time
    return elapsed_initial_time, elapsed_iterations_time


# alleles_count_ = 100
# population_amount_ = 10000000
# interval_for_alpha = 0.125
#
# alpha_values = np.arange(start=0.0, stop=1.0 + interval_for_alpha, step=interval_for_alpha)
# plot_means = []
# plot_stds = []
#
# index = [f'{alpha}' for alpha in alpha_values]
# columns = ['mean', 'std']
# df_positives_for_alpha_values = pd.DataFrame(columns=columns, index=index, dtype=float)
#
# for i, alpha in enumerate(alpha_values):
#     values = np.zeros(20)
#     for experiment in range(20):
#         val = perform_experiment(alleles_count_, population_amount_, alpha, experiment_num=experiment)
#         values[experiment] = val
#     plot_means.append(np.mean(values))
#     plot_stds.append(np.std(values))
#     df_positives_for_alpha_values.iloc[i, 0] = np.mean(values)
#     df_positives_for_alpha_values.iloc[i, 1] = np.std(values)
#
# df_positives_for_alpha_values.to_csv(f'data/positives_for_alpha_values_alleles={alleles_count_}_population={population_amount_}')
# plt.errorbar(alpha_values, plot_means, plot_stds, fmt='o')
# plt.xlabel('Alpha values')
# plt.ylabel('% Positives')
# plt.title(f'Alleles: {alleles_count_}. Population: {population_amount_}')
# plt.show()



if __name__ == '__main__':
    alleles_count_ = 10
    # alleles_count_ = [10, 50, 100, 500, 1000]
    population_amount_ = [1000000, 10000000, 100000000]
    interval_for_alpha = 0.125

    alpha_values = np.arange(start=0.0, stop=1.0 + interval_for_alpha, step=interval_for_alpha)
    plot_means = []
    plot_stds = []

    index = [f'{alpha}' for alpha in alpha_values]
    columns = population_amount_
    df_population_alphas_means = pd.DataFrame(columns=columns, index=index, dtype=float)

    for i, alpha in enumerate(alpha_values):
        for j, population in enumerate(population_amount_):
            values = np.zeros(20)
            for experiment in range(20):
                val = perform_experiment(alleles_count_, population, alpha, experiment_num=experiment)
                values[experiment] = val
            df_population_alphas_means.iloc[i, j] = np.mean(values)

    # df_population_alphas_means.to_csv(f'data/populations_alphas_means_alleles={alleles_count_}')
    sns.heatmap(df_population_alphas_means, annot=True, cmap=plt.cm.CMRmap_r)
    plt.xlabel('Population size')
    plt.ylabel('Alpha values')
    plt.title(f'Means of % positive results. {alleles_count_} Alleles.')
    plt.show()







# alleles_count_ = 100
# population_amount_ = 1000000
# alpha_val_ = 0.0
#
# # if we put 10000 alleles and 100000 population (10 ratio) is problematic with infinite loops
# # number of alleles x amount of population
# alleles_nums = np.array([10, 50, 100, 1000, 5000])
# population_nums = np.array([1000000, 10000000, 100000000])
# experiments_num = 20
# # alleles_nums = np.array([100, 50])
# # population_nums = np.array([1000000, 10000000])
# # experiments_num = 20
#
# index = [f'{alleles}' for alleles in alleles_nums]
# columns = [f'{population}' for population in population_nums]
# df_means = pd.DataFrame(columns=columns, index=index, dtype=float)
# df_stds = pd.DataFrame(columns=columns, index=index, dtype=float)
#
# for i, alleles in enumerate(alleles_nums):
#     for j, population in enumerate(population_nums):
#         experiments_values = np.zeros(experiments_num)
#         for experiment in range(experiments_num):
#             experiments_values[experiment] = perform_experiment(alleles_count=alleles,
#                                                                 population_amount=population,
#                                                                 alpha_val=alpha_val_,
#                                                                 experiment_num=experiment)
#         df_means.iloc[i, j] = np.mean(experiments_values)
#         df_stds.iloc[i, j] = np.std(experiments_values)
#
# # print(initial_times_df)
# # print(iterations_times_df)
# df_means.to_csv(f'data/means_alpha=00')
# df_stds.to_csv(f'data/stds_alpha=00')
