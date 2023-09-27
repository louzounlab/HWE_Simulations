import numpy as np
import random

import pandas as pd
import seaborn as sns
import array as arr

import utils_with_certainty
import matplotlib.pyplot as plt
import time


# populations_probs of size alleles_count x alleles_count x population_amount
# returns: calculated N, {p(i)}, {p(i,j)}, {O(i,j)}
def prepare_experiment_data(alleles_count, population_amount, alpha_val):
    print(f'alleles: {alleles_count}, population: {population_amount}')
    start_time = time.time()

    # probabilities {p(i)}
    alleles_probabilities = utils_with_certainty.calculate_alleles_probabilities(alleles_count)

    # probabilities: {p(i,j)}. only upper triangle
    probabilities = utils_with_certainty.calc_p_ij(alleles_count, alleles_probabilities)
    # print(probabilities)

    # observed: {O(i,j) = Poisson()}. only upper triangle. here alpha affects the non-diagonal elements
    observed = utils_with_certainty.calc_O_ij(population_amount, alleles_count, alpha_val, alleles_probabilities,
                                              probabilities)
    print('calculated observed')

    # calculate new amount of population based on the observations
    population_amount_calculated = 0
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            population_amount_calculated += observed[i, j]
    end_time = time.time()
    elapsed_initial_time = end_time - start_time
    return population_amount_calculated, alleles_probabilities, probabilities, observed, elapsed_initial_time


# given the simulated data (or real data), perform Gibbs Sampling.
# returns:  result (1 for significance or 0)
def perform_experiment(alleles_count, population_amount, alpha_val, experiment_num,
                       population_amount_calculated, alleles_probabilities, probabilities, observed,
                       plot_index=0):
    # [ln(p_0), ln(p_1),...,]
    # actually [0, delta_1, delta_2,...,delta_k]
    list_probabilities = []

    # print(observed)

    # calculate cdf. [(p_11, 1), (p_11+p_12, 2), ..., (p_11+...+p_kk, k)]
    cdf_dict = {}
    utils_with_certainty.calculate_cdf_dict(alleles_count, observed, cdf_dict=cdf_dict)

    # print(cdf_dict)

    # add first ln(probability) as 0. We only care about the deltas anyway.
    list_probabilities.append(0)

    # calc couples
    # first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
    # first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
    i = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                       k=1)[0]
    j = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
                                                      left_alleles=[i, i])
    first_couple = [i, j]

    k = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                       k=1)[0]
    # while k in {i, j}:
    #     k = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                        k=1)[0]
    l = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                       k=1)[0]
    # l = k
    # while l in {i, j}:
    #     l = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                        k=1)[0]
    second_couple = [k, l]
    # pick two alleles using the cdf (first we pick a number between 0 and 1 and then get element with the closest
    # probability)
    # t
    allele_1 = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                              k=1)[0]
    # while allele_1 in {first_couple[0], first_couple[1], second_couple[1]}:
    #     allele_1 = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                               k=1)[0]
    # m
    allele_2 = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
                                                             left_alleles=second_couple)
    # while allele_2 in {allele_1}:
    #     allele_2 = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
    #                                                              left_alleles=second_couple)

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
                                                          population_amount_calculated,
                                                          probabilities)
    # now keep iterating to fill the list of deltas
    # iterations: 100000
    iterations = 5000000 * 10
    # iterations = 5
    for k in range(iterations):
        # sums.append(sum_loops[0])
        if k % 1000 == 0:
            print(f' loop {k} / {iterations}. experiment num: {experiment_num}. alleles: {alleles_count}. '
                  f'population: {population_amount}. alpha: {alpha_val}')
            # print(observed)
            # print('cdf dict:')
            # print(cdf_dict)
        # print(alleles_state)

        # updating the couples
        first_couple = [first_couple[0], couple_from_cdf[1]]
        second_couple = [second_couple[0], couple_from_cdf[0]]

        # pick two alleles using the cdf (first we pick a number between 0 and 1 and then get element with the closest
        # probability)
        # x
        x = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                           k=1)[0]
        # while x in {first_couple[0], first_couple[1], second_couple[1]}:
        #     x = random.choices(population=range(alleles_count), weights=alleles_probabilities,
        #                        k=1)[0]
        # y
        y = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
                                                          left_alleles=second_couple)
        # while y in {x}:
        #     y = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
        #                                                       left_alleles=second_couple)

        couple_from_cdf = [x, y]

        couples = np.zeros(shape=(2, 3), dtype=int)
        couples[0, :] = [first_couple[1], x, +1]
        couples[1, :] = [second_couple[1], y, -1]

        # print(f'plus couple to update: {first_couple[1]}, {x}')
        # print(f'minus couple to update: {second_couple[1]}, {y}')

        utils_with_certainty.update_current_delta_probability(list_probabilities, observed, alleles_probabilities,
                                                              couples,
                                                              population_amount_calculated,
                                                              probabilities)
        # print(observed)

    # now we have the list of probabilities. check if 95% of the elements (sum of deltas) are bigger than 1.
    sum_current = 0
    bigger_counter = 0
    # start_from = 10000
    start_from = utils_with_certainty.calculate_start_time(alleles_count=alleles_count,
                                                           population_amount=population_amount,
                                                           alleles_probabilities=alleles_probabilities,
                                                           observed=observed,
                                                           cdf_dict=cdf_dict)

    values = []

    for delta in list_probabilities[:start_from]:
        sum_current += delta

    for delta in list_probabilities[start_from:]:
        sum_current += delta
        values.append(sum_current)
        if sum_current >= list_probabilities[0]:
            bigger_counter += 1
    # if (bigger_counter / len(list_probabilities) - 1) > 0.95:
    #     return 1

    # print(observed)
    # plt.subplot(7, 7, plot_index)
    print(f'start time: {start_from}, original is 10000, amount of iterations: {iterations}')
    plt.plot(list(range(start_from, len(values) + start_from)), values)
    plt.show()
    # plt.xlabel('Time')
    # plt.ylabel('Approximated ln(p_i)')
    # plt.title(f'Alleles: {alleles_count}. Population: {population_amount}')

    result = bigger_counter / (len(list_probabilities) - start_from)
    print(result)
    return result


###### TEST
# (alleles_count, population_amount, alpha_val, experiment_num,
#                        population_amount_calculated, alleles_probabilities, probabilities, observed)
plt.rcParams["font.family"] = "Arial"

# alleles_count_ = 4
# population_amount_ = 1000000
# alpha_val_ = 0.0
# experiment_num_ = 0
# #alleles_probabilities_ = np.zeros(alleles_count_) + 0.25
#
# # np.seterr(all='raise')
# plt.subplot (7, 7, 1)
# for i in range(1, 50):
#     population_amount_calculated_, alleles_probabilities, probabilities_, observed_, elapsed_initial_time = \
#         prepare_experiment_data(alleles_count=alleles_count_,
#                                 population_amount=population_amount_,
#                                 alpha_val=alpha_val_)
#
#     print(f'population amount calculated = {population_amount_calculated_}')
#     perform_experiment(alleles_count_, population_amount_, alpha_val_, experiment_num_,
#                        population_amount_calculated_, alleles_probabilities, probabilities_, observed_, i)
# plt.show()

alleles_count_ = 4
population_amount_ = 1000000
alpha_val_ = 0.0
experiment_num_ = 0

population_amount_calculated_, alleles_probabilities, probabilities_, observed_, elapsed_initial_time = \
    prepare_experiment_data(alleles_count=alleles_count_,
                            population_amount=population_amount_,
                            alpha_val=alpha_val_)

perform_experiment(alleles_count_, population_amount_, alpha_val_, experiment_num_,
                   population_amount_calculated_, alleles_probabilities, probabilities_, observed_)
