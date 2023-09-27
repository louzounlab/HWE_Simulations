import random
import numpy as np
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
def calculate_marginal_probabilities(alleles_count):
    probs = np.zeros(shape=(alleles_count, alleles_count))
    for i_ in range(alleles_count):
        marginal = calculate_alleles_probabilities(alleles_count)
        probs[:, i_] = marginal
    return probs


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


def calculate_probability_given_i_j(marginal_probabilities_, i_, j_, k_1, k_2):
    prob = marginal_probabilities_[k_1, i_] * marginal_probabilities_[k_2, j_]
    # EVEN
    if i_ != j_:
        prob += marginal_probabilities_[k_1, j_] * marginal_probabilities_[k_2, i_]
    return prob


# calculate the denominator for the new Chi-Squared test
def calculate_total_variance_alleles_i_j(alleles_amount_,
                                         population_amount,
                                         alleles_probabilities, marginal_probabilities_,
                                         uncertainty_,
                                         observed_,
                                         sum_observed_,
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
    expected_probability = alleles_probabilities[i_] * alleles_probabilities[j_] * mult
    expected = population_amount * alleles_probabilities[i_] * alleles_probabilities[j_] * mult

    i_j_observed_probability = observed_[i_, j_] / sum_observed_

    return expected + expected_probability * uncertainty_ * (population_amount ** 1.5)


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


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   should_use_new_test_):
    # matrix where every row is the alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, alleles_count), dtype=np.int8)

    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)

    # probabilities: A(j,i) = p_i(j), columns i
    marginal_probabilities = calculate_marginal_probabilities(alleles_count)

    # print(f'alleles: {alleles_probabilities}')
    # print(f' marginal: {marginal_probabilities}')

    # 1) assignment of alleles with certainty
    for k in range(population_amount):
        # choice allele according to alpha, if choice is zero then regular like model A
        choice = random.choices(population=[0, 1], weights=[alpha_val, 1 - alpha_val], k=1)[0]
        # print(f'choice: {choice}, k={k}')
        # alleles = np.zeros(shape=(2, 0))
        if choice == 0:
            alleles = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                                     k=alleles_count)
        else:
            allele = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                                    k=1)[0]
            alleles = [allele] * alleles_count
        alleles = np.array(alleles)
        alleles_individuals[k, :] = alleles

        # print(f'alleles_individuals: {alleles_individuals}')

    # adding uncertainty to our model
    # DEBUG
    for k in range(population_amount):  # for each person k
        choice = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=1)[0]
        if choice == 0:
            for i in range(alleles_count):  # represent an index for the set of alleles
                allele_i = alleles_individuals[k, i]  # the i-th allele of person k
                # the new allele we will choose from the marginal probability of allele_i
                new_allele = \
                    random.choices(population=range(alleles_count), weights=marginal_probabilities[:, allele_i],
                                   k=1)[0]
                # assign the new allele for the person k
                alleles_individuals[k, i] = new_allele

    # print(f'alleles_individuals after uncertainty: {alleles_individuals}')

    # count the amounts of alleles occurrences in the population
    counts_observed = np.zeros((alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # count amount of occurrences of j,k
            count_tuple = count_row_occurrence_in_2d_array(j, k, alleles_individuals)

            counts_observed[j, k] = count_tuple
            counts_observed[k, j] = count_tuple

    counts_expected = np.zeros((alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # EVEN
            mult = 1
            if k != j:
                mult = 2
            expected_value = mult * population_amount * alleles_probabilities[k] * alleles_probabilities[j]
            counts_expected[k, j] = expected_value
            counts_expected[j ,k] = expected_value

    sum_observed = 0.0
    for k in range(alleles_count):
        for j in range(k, alleles_count):
            sum_observed += counts_observed[k, j]

    counts_total_variance = np.zeros((alleles_count, alleles_count))
    for k in range(alleles_count):
        for j in range(k, alleles_count):
            total_variance = calculate_total_variance_alleles_i_j(alleles_count,
                                                                  population_amount,
                                                                  alleles_probabilities, marginal_probabilities,
                                                                  uncertainty_val,
                                                                  counts_observed,
                                                                  sum_observed,
                                                                  k, j)

            counts_total_variance[k, j] = total_variance
            counts_total_variance[j, k] = total_variance

    chi_squared_stat = calculate_chi_squared_value(counts_expected_=counts_expected,
                                                   counts_observed_=counts_observed,
                                                   counts_total_variance_=counts_total_variance,
                                                   should_use_new_test_=should_use_new_test_)
    dof = (alleles_count * (alleles_count + 1)) / 2 - 1

    # print(f' alpha for choice: {alpha_val}')
    # print(f' chi square value: {chi_squared_stat}')

    crit = stats.chi2.ppf(q=0.95, df=dof)
    # print(f'Critical value: {crit}')

    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,
                                 df=dof)
    if p_value < 0.05:
        return 1
    return 0


fix, ax = plt.subplots()
ax.grid()
ax.set_axisbelow(True)

num_of_experiments = 400  # 500 amount of experiments for each alpha ranging from zero to one
alleles_amount = 2  # 2
population_size = 20000  # 35
interval_for_alpha = 0.04  # 0.02
# uncertainty = 0.2
uncertainty_vals = [0.0, 0.1, 0.2, 0.4]
# blue color: 'royalblue'
colors = ['royalblue', 'slategrey', 'limegreen', 'deeppink']
#alpha_values = np.arange(start=0, stop=1 + interval_for_alpha, step=interval_for_alpha)  # start, stop, step
alpha_values = np.array([1.0])
confidence_amounts_per_alpha = np.zeros((len(uncertainty_vals), alpha_values.shape[0]))
should_use_new_test = 0
# markers = ['', '>']
markers = ['.', '>']
labels = ['Old Chi-Squared', 'Corrected Chi-Squared']
# print(alpha_values.shape)
# print(confidence_amounts_per_alpha.shape)
for should_use_new_test in range(2):
    for uncertainty_index, uncertainty in enumerate(uncertainty_vals):
        for i, alpha in enumerate(alpha_values):
            counter = 0

            for experiment in range(num_of_experiments):
                counter += run_experiment(alleles_amount, population_size, alpha, uncertainty,
                                          should_use_new_test)
            print(f'uncertainty = {uncertainty}, alpha = {alpha}, test type = {should_use_new_test}')
            confidence_amount = counter / num_of_experiments
            confidence_amounts_per_alpha[uncertainty_index, i] = confidence_amount

        print(confidence_amounts_per_alpha[uncertainty_index, :])
        plt.plot(alpha_values + should_use_new_test, confidence_amounts_per_alpha[uncertainty_index, :],
                 label=f'{labels[should_use_new_test]}, uncertainty = {uncertainty}',
                 marker=markers[should_use_new_test],
                 color=colors[uncertainty_index])

# print(confidence_amounts_per_alpha[-1])

# plt.plot(alpha_values, confidence_amounts_per_alpha, 'r-', label="")
plt.xlabel('Alpha values')
plt.ylabel('Percentage of confidence amounts')
plt.title('Old vs Corrected Chi-Squared')
plt.legend()
plt.savefig('chi_squared_simulation', pad_inches=0.2, bbox_inches="tight")
plt.show()