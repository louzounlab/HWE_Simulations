import random
import numpy as np
import seaborn as sns
import math
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt


#####
# Special comments:
# EVEN: means we take care of the case i!=j different from the case i=j
# DEBUG: didn't debug this part yet
#####


# assuming x is 2-dimensional
def softmax_1d(x):
    sum = 0.0
    for row in range(x.shape[0]):
        sum += np.exp(x[row])
    return np.exp(x) / sum


# generate a vector of probabilities
def calculate_alleles_probabilities(alleles_count):
    probs = np.random.uniform(0.0, 2.0, size=alleles_count)
    probs = softmax_1d(probs)
    return probs


# generate a marginal probabilities matrix
# column i represent the probabilities p_i(j), j are for the rows
# def calculate_marginal_probabilities(alleles_count):
#     probs = np.zeros(shape=(alleles_count, alleles_count))
#     for i_ in range(alleles_count):
#         marginal = calculate_alleles_probabilities(alleles_count)
#         probs[:, i_] = marginal
#     return probs


# data: N x 2, row = (j, k). alleles j,k
def count_row_occurrence_in_2d_array(j, k, data):
    count = 0
    # for every row (person) in population
    for row in range(data.shape[0]):
        if (data[row, 0] == j and data[row, 1] == k) or (data[row, 0] == k and data[row, 1] == j):
            count += 1
            # count += (alleles_probabilities[j] * alleles_probabilities[k] / (population_probs[j,k,row]))
            # count += 1

    return count


# def calculate_probability_given_i_j(marginal_probabilities_, i_, j_, k_1, k_2):
#     prob = marginal_probabilities_[k_1, i_] * marginal_probabilities_[k_2, j_]
#     # EVEN
#     if i_ != j_:
#         prob += marginal_probabilities_[k_1, j_] * marginal_probabilities_[k_2, i_]
#     return prob


# calculate the denominator for the new Chi-Squared test
def calculate_total_variance_alleles_i_j(alleles_amount_,
                                         population_amount,
                                         alleles_probabilities,
                                         observed_,
                                         i_, j_):
    # mult = 1
    # if i_ != j_:
    #     mult = 2
    # expected = population_amount * alleles_probabilities[i_] * alleles_probabilities[j_] * mult
    #
    # k, l = min(i_, j_), max(i_, j_)
    # p_kl = observed_[k, l] / ((alleles_amount_ * (alleles_amount_ - 1)) / 2)
    # U = 0
    # for allele in range(alleles_amount_):
    #     U -= (alleles_probabilities[allele] ** 2)
    # return expected + p_kl * U

    # EVEN
    # mult = 1
    # if i_ != j_:
    #     mult = 2
    # expected = population_amount * alleles_probabilities[i_] * alleles_probabilities[j_] * mult
    # probability_not_i_j = 0
    # for k_1 in range(alleles_amount):
    #     for k_2 in range(k_1, alleles_amount):
    #         if {k_1, k_2} != {i_, j_}:
    #             probability_not_i_j += marginal_probabilities_[k_1, i_] * marginal_probabilities_[k_2, j_]
    #             # EVEN
    #             if i_ != j_:
    #                 probability_not_i_j += marginal_probabilities_[k_1, j_] * marginal_probabilities_[k_2, i_]
    #
    # probability_i_j = 0
    # denominator = 0
    # for k_1 in range(alleles_amount_):
    #     for k_2 in range(k_1, alleles_amount_):
    #         if {k_1, k_2} != {i_, j_}:
    #             probability_i_j += calculate_probability_given_i_j(marginal_probabilities_, k_1, k_2,
    #                                                                i_, j_)
    #             denominator += alleles_probabilities[k_1] * alleles_probabilities[k_2]
    #             # EVEN
    #             if k_1 != k_2:
    #                 denominator += alleles_probabilities[k_1] * alleles_probabilities[k_2]
    #
    # probability_i_j = probability_i_j / denominator
    #
    # # probability_not_i_j = 1 - marginal_probabilities_[i_, i_] * marginal_probabilities_[j_, j_]
    # # if i_ != j_:
    # #     probability_not_i_j -= marginal_probabilities_[j_, i_] * marginal_probabilities_[i_, j_]
    # var = expected + expected * uncertainty_ * 1 * probability_not_i_j + \
    #     population_amount * (1 - alleles_probabilities[i_] * alleles_probabilities[j_] * mult) * uncertainty_ * 0.3 * probability_i_j

    # EVEN
    mult = 1
    if i_ != j_:
        mult = 2
    # expected_probability = alleles_probabilities[i_] * alleles_probabilities[j_] * mult
    expected = population_amount * alleles_probabilities[i_] * alleles_probabilities[j_] * mult

    # i_j_observed_probability = observed_[i_, j_] / sum_observed_

    return expected


# need to calculate only on the upper triangle because the matrices are symmetric
def calculate_chi_squared_value(counts_expected_, counts_observed_, counts_total_variance_, should_use_new_test_):
    value = 0
    for row in range(counts_expected_.shape[0]):
        for col in range(row, counts_expected_.shape[1]):
            expected_ = counts_expected_[row, col]
            observed_ = counts_observed_[row, col]
            variance_ = counts_total_variance_[row, col]
            if should_use_new_test_:
                denominator = variance_
            else:
                denominator = expected_
            value += (((expected_ - observed_) ** 2) / denominator)
    return value


# A i_hat,j_hat = p_i(i_hat) * p_j(j_hat). each element (k,l) represent to probability to get (k,l) alleles
# given I had (i,j) before
# def calculate_uncertainty_matrix(alleles_count, marginal_probabilities, i, j):
#     uncertainty_matrix = np.zeros(shape=(alleles_count, alleles_count))
#     for i_hat in range(alleles_count):
#         for j_hat in range(alleles_count):
#             uncertainty_matrix[i_hat, j_hat] = marginal_probabilities[i_hat, i] * marginal_probabilities[j_hat, j]
#             if i != j:  # EVEN
#                 # we add the probability that we moved from (j,i) to (i_hat, j_hat)
#                 uncertainty_matrix[i_hat, j_hat] += marginal_probabilities[i_hat, j] * marginal_probabilities[j_hat, i]
#
#     return uncertainty_matrix


def calc_variance(input_list, input_mean):
    var = 0.0
    for element in input_list:
        var += (element - input_mean) ** 2
    return var


def prepare_probabilities(alleles_count):
    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)

    return alleles_probabilities


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   alleles_probabilities, plot_index):
    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (
                        alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
                probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]

    # matrix where every row is the alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, 2), dtype=np.int32)

    # print(f'alleles: {alleles_probabilities}')
    # print(f' marginal: {marginal_probabilities}')

    # 1) assignment of alleles with certainty
    probabilities_list = probabilities.flatten()
    indices = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=population_amount)
    for k in range(population_amount):
        # notice the alleles i != j have probability 2 * p(i,j)
        index = indices[k]
        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

    # now we multiply the upper triangle by 2
    # for t in range(alleles_count):
    #     for m in range(t + 1, alleles_count):
    #         probabilities[t, m] *= 2
    #         probabilities[m, t] = 0

    # print(f'alleles_individuals: {alleles_individuals}')

    # matrix p_k(i,j) = A[i,j,k]
    all_probabilities = np.zeros(shape=(alleles_count, alleles_count, population_amount))

    choices = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=population_amount)

    # adding uncertainty to our model
    for k in range(population_amount):
        # person k has alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)

        # choice whether this person will have uncertain alleles
        choice = choices[k]

        # this person has certain alleles
        if choice == 1:
            all_probabilities[j, l, k] = 1.0
        # this person has uncertain alleles
        if choice == 0:
            all_probabilities[:, :, k] = probabilities

    # now we have the probabilities {p_k(i,j)}
    # lets get the final p(i,j) = 1 / N * sum_k p_k(i,j)
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            probability = 0.0
            for k in range(population_amount):
                if t == m:
                    probability += all_probabilities[t, m, k]
                else:
                    probability += (all_probabilities[t, m, k] + all_probabilities[m, t, k])
            probability /= population_amount

            observed_probabilities[t, m] = probability

    observations = observed_probabilities * population_amount

    # sum_observed = 0.0
    # for k in range(alleles_count):
    #     for j in range(k, alleles_count):
    #         sum_observed += observations[k, j]

    x_values = []
    counts_expected = np.zeros((alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # EVEN
            mult = 1
            if k != j:
                mult = 2
            expected_value = mult * population_amount * alleles_probabilities[j] * alleles_probabilities[k]
            counts_expected[j, k] = expected_value
            x_values.append(counts_expected[j, k])
    # counts_total_variance = np.zeros((alleles_count, alleles_count))
    # for j in range(alleles_count):
    #     for k in range(j, alleles_count):
    #         total_variance = calculate_total_variance_alleles_i_j(alleles_count,
    #                                                               population_amount,
    #                                                               alleles_probabilities,
    #                                                               observations,
    #                                                               j, k)
    #
    #         counts_total_variance[j, k] = total_variance

    observed_no_sampling_vector = []
    # variance_vector = []

    for j in range(alleles_count):
        for k in range(j, alleles_count):
            observed_no_sampling_vector.append(observations[j, k])
            # variance_vector.append(counts_total_variance[j, k])

    # observed_vector = [np.log(np.log(observed_value)) for observed_value in observed_vector]
    # observed_vector = [np.log(np.log(observed_value)) for observed_value in observed_vector]
    #
    # observed_diagonal_vector = [np.log(np.log(observed_value)) for observed_value in observed_diagonal_vector]
    # variance_diagonal_vector = [np.log(np.log(observed_value)) for observed_value in variance_diagonal_vector]
    plt.subplot(2, 3, plot_index)
    plt.scatter(x_values, observed_no_sampling_vector, color='black', label='OBS NO SAMPLING', s=1.5)
    plt.xlabel('N*p(i)*p(j)')

    # for each k make {p_k(i,j)} symmetric and choose for every k alleles i,j
    for k in range(population_amount):
        # for t in range(alleles_count):
        #     for m in range(t, alleles_count):
        #         all_probabilities[t, m, k] = (all_probabilities[t, m, k] + all_probabilities[m, t, k]) / 2
        probabilities_list = all_probabilities[:, :, k].flatten()
        index = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=1)[0]

        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

    observed_sampling_vector = []
    variance_sampling_vector = []
    variance_no_sampling_vector = []

    # calculate observations
    observations = np.zeros(shape=(alleles_count, alleles_count))
    for k in range(population_amount):
        t, m = alleles_individuals[k]
        observations[t, m] += 1

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            observed_sampling_vector.append(observations[t, m])
    # calculate variance

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            probabilities_list = all_probabilities[t, m, :].flatten()
            mean_value = np.mean(probabilities_list)
            variance = calc_variance(input_list=probabilities_list,
                                     input_mean=mean_value)
            variance_no_sampling_vector.append(variance)

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            indicators = []
            for k in range(population_amount):
                y, u = alleles_individuals[k]
                indicator = int((t == y) and (m == u))
                indicators.append(indicator)
            mean_value = np.mean(indicators)
            variance = calc_variance(input_list=indicators,
                                     input_mean=mean_value)
            variance_sampling_vector.append(variance)

    # observed_vector = [np.log(np.log(observed_value)) for observed_value in observed_vector]
    # observed_vector = [np.log(np.log(observed_value)) for observed_value in observed_vector]
    #
    # observed_diagonal_vector = [np.log(np.log(observed_value)) for observed_value in observed_diagonal_vector]
    # variance_diagonal_vector = [np.log(np.log(observed_value)) for observed_value in variance_diagonal_vector]
    # min_x2 = min(observed_vector + observed_diagonal_vector)
    # max_x2 = max(observed_vector + observed_diagonal_vector)
    #
    # min_y2 = min(variance_vector + variance_diagonal_vector)
    # max_y2 = max(variance_vector + variance_diagonal_vector)
    #
    # min_x = min(min_x, min_x2)
    # max_x = max(max_x, max_x2)
    #
    # min_y = min(min_y, min_y2)
    # max_y = max(max_y, max_y2)
    #
    # min_val = min(min_x, min_y)
    # max_val = max(max_x, max_y)
    plt.scatter(x_values, observed_sampling_vector, color='orange', label='OBS SAMPLING', s=1.5)
    plt.scatter(x_values, variance_no_sampling_vector, color='blue', label='VAR NO SAMPLING', s=1.5)
    plt.scatter(x_values, variance_sampling_vector, color='lime', label='VAR SAMPLING', s=1.5)
    # plt.xlim([min_val, max_val])
    # plt.ylim([min_val, max_val])
    plt.legend()
    plt.title(f'Alpha: {alpha_val}. Uncertainty: {uncertainty_val}')


if __name__ == '__main__':
    alleles_amount = 50  # 20
    population_size = 10000  # 10000
    alpha_vals = [1.0, 0.85]
    uncertainty_vals = [0.0, 0.2, 0.4]

    sns.set_style('white')
    plt.rcParams["font.family"] = "Arial"

    fig, axes = plt.subplots(len(alpha_vals), len(uncertainty_vals), constrained_layout=True)

    print('Calculating probabilities')
    # get {p(i)}, {p(i|j)}
    alleles_probabilities_ = prepare_probabilities(alleles_count=alleles_amount)

    plot_num = 0

    for alpha in alpha_vals:
        for uncertainty in uncertainty_vals:
            print(f'alpha: {alpha}, uncertainty: {uncertainty}')
            plot_num += 1
            run_experiment(alleles_count=alleles_amount,
                           population_amount=population_size,
                           alpha_val=alpha, uncertainty_val=uncertainty,
                           alleles_probabilities=alleles_probabilities_,
                           plot_index=plot_num)
    plt.savefig('chi_squared.png', pad_inches=0.2, bbox_inches="tight")
    plt.show()
