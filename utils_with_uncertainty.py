import numpy as np


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
                                         alleles_probabilities,
                                         uncertainty_,
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
    expected = population_amount * alleles_probabilities[i_] * alleles_probabilities[j_] * mult

    # i_j_observed_probability = observed_[i_, j_] / sum_observed_
    i_j_observed_probability = observed_[i_, j_]

    return expected + i_j_observed_probability * uncertainty_ # * (population_amount ** 1.5)


# need to calculate only on the upper triangle because the matrices are symmetric
def calculate_chi_squared_value(alleles_amount, population_amount, alleles_probabilities, counts_observed_, counts_total_variance_, should_use_new_test_):
    value = 0
    for row in range(alleles_amount):
        for col in range(alleles_amount):
            mult = 1
            if row != col:
                mult = 2
            expected = population_amount * alleles_probabilities[row] * alleles_probabilities[col] * mult
            observed_ = counts_observed_[row, col]
            variance_ = counts_total_variance_[row, col]
            if should_use_new_test_:
                denominator = variance_
            else:
                denominator = expected
            value += (((expected - observed_) ** 2) / denominator)
    return value
