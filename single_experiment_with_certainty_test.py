import numpy as np
import random

import pandas as pd
import seaborn as sns
import array as arr

import utils_with_certainty
import matplotlib.pyplot as plt
import time
import warnings


def get_O_ij(observed, i, j):
    row, col = min(i, j), max(i, j)
    return observed[row, col]


def get_p_ij(probabilities, i, j):
    i, j = min(i, j), max(i, j)
    return probabilities[i, j]


def calc_observed_cdf(alleles_count, observed):
    # list of numpy arrays
    observed_cdf = []
    for i in range(alleles_count):
        np_array = np.zeros(alleles_count)
        # 1,i 2,i 3,i ,..., i,i i,i+1 i,i+2 ,..., i,n
        for k in range(alleles_count):
            t, m = min(i, k), max(i, k)
            mult = 0.5
            if t == m:
                mult = 1.0
            np_array[k] = mult * observed[t, m]
        observed_cdf.append(np_array)
    return observed_cdf


def update_current_delta_probability(list_probabilities: list, observed, observed_cdf, alleles_probabilities, couples,
                                     population_amount_calculated, probabilities):
    for row in range(couples.shape[0]):
        # here we make sure that i <= j so we can access the upper triangle of probabilities
        i, j = min(couples[row, 0], couples[row, 1]), max(couples[row, 0], couples[row, 1])
        # i, j = int(i), int(j)

        # z = np.log(population_amount_calculated * probabilities[i, j] / (1 - probabilities[i, j]))
        sign = couples[row, 2]
        val = 0
        try:
            if sign == -1:
                # val = np.log(observed[i, j] / (population_amount_calculated - observed[i, j] + 1)) \
                #       + np.log((1 - get_p_ij(alleles_probabilities, i, j)) / get_p_ij(alleles_probabilities, i, j))
                val = np.log(get_O_ij(observed, i, j)) - np.log(
                    population_amount_calculated - get_O_ij(observed, i, j) + 1) \
                      - np.log(get_p_ij(probabilities, i, j)) + np.log(1 - get_p_ij(probabilities, i, j))
            elif sign == 1:
                val = np.log(population_amount_calculated - get_O_ij(observed, i, j)) - np.log(
                    get_O_ij(observed, i, j) + 1) \
                      + np.log(get_p_ij(probabilities, i, j)) - np.log(1 - get_p_ij(probabilities, i, j))
        except RuntimeWarning:
            print(f'''
i,j: {i}, {j}.
O_ij: {get_O_ij(observed, i, j)}
population calculated: {population_amount_calculated}
p_ij: {get_p_ij(probabilities, i, j)}
sign: {sign}
observations: {observed}
alleles probabilities: {alleles_probabilities}''')
            raise ValueError('A very specific bad thing happened.')

        # print(observed)
        # print(f'i: {i}, j: {j}')
        # sum_current += sign * (z + np.log(observed[i, j]))
        list_probabilities.append(val)
        # unbalanced_alleles[i] += sign
        # unbalanced_alleles[j] += sign
        observed[i, j] += sign

        # also update observed_cdf
        # o(i,j)+1 -> j_th row, i_th row
        if i == j:
            observed_cdf[j][i] += sign
        else:
            observed_cdf[j][i] += 0.5 * sign
            observed_cdf[i][j] += 0.5 * sign


# populations_probs of size alleles_count x alleles_count x population_amount
# returns: calculated N, {p(i)}, {p(i,j)}, {O(i,j)}
def prepare_experiment_data(alleles_count, population_amount, alpha_val):
    start_time = time.time()

    # probabilities {p(i)}
    alleles_probabilities = utils_with_certainty.calculate_alleles_probabilities(alleles_count)

    # probabilities: {p(i,j)}. only upper triangle
    probabilities = utils_with_certainty.calc_p_ij(alleles_count, alleles_probabilities)
    # print(probabilities)

    # observed: {O(i,j) = Poisson()}. only upper triangle. here alpha affects the non-diagonal elements
    observed = utils_with_certainty.calc_O_ij(population_amount, alleles_count, alpha_val, alleles_probabilities,
                                              probabilities)
    observed_cdf = calc_observed_cdf(alleles_count, observed)
    print('calculated observed')

    # calculate new amount of population based on the observations
    population_amount_calculated = 0
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            population_amount_calculated += observed[i, j]
    end_time = time.time()
    elapsed_initial_time = end_time - start_time
    return population_amount_calculated, alleles_probabilities, probabilities, observed, observed_cdf, elapsed_initial_time


# given the simulated data (or real data), perform Gibbs Sampling.
# returns:  result (1 for significance or 0)
def perform_experiment(alleles_count,
                       population_amount_calculated, alleles_probabilities, probabilities, observed, observed_cdf,
                       plot_index=0):
    # [ln(p_0), ln(p_1),...,]
    # actually [0, delta_1, delta_2,...,delta_k]
    list_probabilities = []

    range_alleles_count = range(alleles_count)

    # print(observed)

    # calculate cdf. [(p_11, 1), (p_11+p_12, 2), ..., (p_11+...+p_kk, k)]
    # cdf_dict = {}
    # utils_with_certainty.calculate_cdf_dict(alleles_count, observed, cdf_dict=cdf_dict)

    # print(cdf_dict)

    # add first ln(probability) as 0. We only care about the deltas anyway.
    list_probabilities.append(0)

    # calc couples
    # first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
    # first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
    i = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                       k=1)[0]

    j = random.choices(population=range_alleles_count, weights=observed_cdf[i],
                       k=1)[0]
    first_couple = [i, j]

    k = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                       k=1)[0]
    # while k in {i, j}:
    #     k = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                        k=1)[0]
    l = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                       k=1)[0]
    # l = k
    # while l in {i, j}:
    #     l = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                        k=1)[0]
    second_couple = [k, l]

    # every row represents a couple i,j and a value (either 1 or -1)
    # in the first iteration here we have 4 couples, but in the next iterations we will have 2 couples.
    couples = np.zeros(shape=(2, 3), dtype=int)
    couples[0, :] = [first_couple[0], first_couple[1], -1]
    couples[1, :] = [second_couple[0], second_couple[1], +1]

    # calculate new delta and modify observed
    update_current_delta_probability(list_probabilities, observed, observed_cdf, alleles_probabilities, couples,
                                     population_amount_calculated,
                                     probabilities)

    # pick two alleles using the cdf (first we pick a number between 0 and 1 and then get element with the closest
    # probability)
    # t
    allele_1 = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                              k=1)[0]
    # while allele_1 in {first_couple[0], first_couple[1], second_couple[1]}:
    #     allele_1 = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                               k=1)[0]
    # m
    allele_2 = random.choices(population=range_alleles_count, weights=observed_cdf[second_couple[1]],
                              k=1)[0]
    # while allele_2 in {allele_1}:
    #     allele_2 = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
    #                                                              left_alleles=second_couple)

    couple_from_cdf = [allele_1, allele_2]

    # every row represents a couple i,j and a value (either 1 or -1)
    # in the first iteration here we have 4 couples, but in the next iterations we will have 2 couples.
    couples[0, :] = [first_couple[1], allele_1, +1]
    couples[1, :] = [second_couple[1], allele_2, -1]

    # calculate new delta and modify observed
    update_current_delta_probability(list_probabilities, observed, observed_cdf, alleles_probabilities, couples,
                                     population_amount_calculated,
                                     probabilities)
    # now keep iterating to fill the list of deltas
    # iterations: 100000
    iterations = 100000
    # iterations = 5
    for k in range(iterations):
        # sums.append(sum_loops[0])
        # if (k % 100) == 0:
        #     print(f' loop {k} / {iterations}. experiment num: {experiment_num}. alleles: {alleles_count}. '
        #           f'population: {population_amount}. alpha: {alpha_val}')
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
        y = random.choices(population=range_alleles_count, weights=observed_cdf[second_couple[1]],
                           k=1)[0]
        # while y in {x}:
        #     y = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
        #                                                       left_alleles=second_couple)

        couple_from_cdf = [x, y]

        couples[0, :] = [first_couple[1], x, +1]
        couples[1, :] = [second_couple[1], y, -1]
        # print(couples)

        # print(f'plus couple to update: {first_couple[1]}, {x}')
        # print(f'minus couple to update: {second_couple[1]}, {y}')

        update_current_delta_probability(list_probabilities, observed, observed_cdf, alleles_probabilities,
                                         couples,
                                         population_amount_calculated,
                                         probabilities)
        # print(observed)

    # now we have the list of probabilities. check if 95% of the elements (sum of deltas) are bigger than 1.
    sum_current = 0
    bigger_counter = 0
    # start_from = 10000
    start_from = 10000

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
    # print(observed)
    # print(observed_cdf)
    # # print(f'start time: {start_from}, original is 10000, amount of iterations: {iterations}')
    # plt.plot(list(range(start_from, len(values) + start_from)), values)
    # plt.show()
    # plt.xlabel('Time')
    # plt.ylabel('Approximated ln(p_i)')
    # plt.title(f'Alleles: {alleles_count}. Population: {population_amount}')

    result = bigger_counter / (len(list_probabilities) - start_from)
    p_value = 1.0 - 2 * abs(max(0.5, result) - 0.5)
    return p_value


###### TEST
# (alleles_count, population_amount, alpha_val, experiment_num,
#                        population_amount_calculated, alleles_probabilities, probabilities, observed)
if __name__ == "__main__":
    plt.rcParams["font.family"] = "Arial"

    alleles_count_ = 1000
    population_amount_ = 1000000
    alpha_val_ = 0.0
    experiment_num_ = 0

    population_amount_calculated_, alleles_probabilities, probabilities_, observed_, observed_cdf_, elapsed_initial_time = \
        prepare_experiment_data(alleles_count=alleles_count_,
                                population_amount=population_amount_,
                                alpha_val=alpha_val_)
    observed_copy = np.copy(observed_)
    observed_cdf_copy = np.copy(observed_cdf_)

    plt.subplot(7, 7, 1)
    for i in range(1, 50):

        print(f'plot {i} / 49')

        perform_experiment(alleles_count_, population_amount_calculated_, alleles_probabilities, probabilities_,
                           observed_copy, observed_cdf_copy, i)
        np.copyto(observed_copy, observed_)
        np.copyto(observed_cdf_copy, observed_cdf_)
    plt.show()

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

    # warnings.filterwarnings("error")
    #
    # alleles_count_ = 1000
    # population_amount_ = 1000000
    # alpha_val_ = 0.0
    # experiment_num_ = 0
    #
    # population_amount_calculated_, alleles_probabilities_, probabilities_, observed_, observed_cdf_, elapsed_initial_time = \
    #     prepare_experiment_data(alleles_count=alleles_count_,
    #                             population_amount=population_amount_,
    #                             alpha_val=alpha_val_)
    #
    # perform_experiment(alleles_count_, population_amount_, alpha_val_, experiment_num_,
    #                    population_amount_calculated_, alleles_probabilities_, probabilities_, observed_, observed_cdf_)

    # alleles_count_ = 2
    # population_amount_ = 983
    # alpha_val_ = 0.0
    # experiment_num_ = 0
    # population_amount_calculated_ = 983
    # alleles_probabilities = np.array([0.5, 0.5])
    # probabilities_ = np.array([[0.25, 0.5],
    #                            [0, 0.25]])
    # observed_ = np.array([[576, 407.],
    #                       [0, 0.]])
    # perform_experiment(alleles_count_, population_amount_, alpha_val_, experiment_num_,
    #                    population_amount_calculated_, alleles_probabilities, probabilities_, observed_)
