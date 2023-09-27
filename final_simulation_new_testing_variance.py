import random
import numpy as np
import math
import statistics
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

    return expected + expected_probability * uncertainty_ * (population_amount)


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


def prepare_probabilities(alleles_count):

    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)

    # probabilities: A(j,i) = p_i(j), columns i
    marginal_probabilities = calculate_marginal_probabilities(alleles_count)

    return alleles_probabilities, marginal_probabilities


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   alleles_probabilities, marginal_probabilities):
    # matrix where every row is the alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, alleles_count), dtype=np.int8)

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

    # counts_expected = np.zeros((alleles_count, alleles_count))
    # for j in range(alleles_count):
    #     for k in range(j, alleles_count):
    #         # EVEN
    #         mult = 1
    #         if k != j:
    #             mult = 2
    #         expected_value = mult * population_amount * alleles_probabilities[k] * alleles_probabilities[j]
    #         counts_expected[k, j] = expected_value
    #         counts_expected[j, k] = expected_value
    #
    # sum_observed = 0.0
    # for k in range(alleles_count):
    #     for j in range(k, alleles_count):
    #         sum_observed += counts_observed[k, j]

    return counts_observed


plt.rcParams["font.family"] = "Arial"
fig, axes = plt.subplots(2,2, constrained_layout=True)

num_of_experiments = 100  # 500 amount of experiments for each alpha ranging from zero to one
alleles_amount = 20  # 20
population_size = 10000  # 10000
alpha_vals = [1.0, 0.8]
uncertainty_vals = [0.0, 0.2]



print('Calculating probabilities')
# get {p(i)}, {p(i|j)}
alleles_probabilities_, marginal_probabilities_ = prepare_probabilities(alleles_count=alleles_amount)
print('Calculating Variance matrix')
variance_matrix = np.zeros(shape=(alleles_amount, alleles_amount))

plot_num = 0

for alpha in alpha_vals:
    for uncertainty in uncertainty_vals:
        plot_num += 1
        # build variance matrix
        for i in range(alleles_amount):
            for j in range(i, alleles_amount):
                mult = 1.0
                if i != j:
                    mult = 2.0
                variance = population_size * alleles_probabilities_[i] * alleles_probabilities_[j] * mult + \
                    alleles_probabilities_[i] * alleles_probabilities_[j] * mult * uncertainty
                variance_matrix[i, j] = variance

        # define a list where every element is a matrix {O(i,j)}
        list_of_observations = []

        print('Running experiments')

        for experiment in range(num_of_experiments):
            print(f'Experiment num: {experiment} / {num_of_experiments}: plot: {plot_num} / {len(alpha_vals) * len(uncertainty_vals)}')
            observations = run_experiment(alleles_count=alleles_amount, population_amount=population_size,
                                          alpha_val=alpha, uncertainty_val=uncertainty, alleles_probabilities=alleles_probabilities_,
                                          marginal_probabilities=marginal_probabilities_)
            list_of_observations.append(observations)

        # lists for the plot
        sample_variances = []
        variances = []

        print('Calculating sample variances for plot')

        # for each i<=j calculate sample variance of {O(i,j)} over the experiments
        for i in range(alleles_amount):
            for j in range(i, alleles_amount):
                # for the specific i,j get all the O(i,j) from all the experiments, as a list
                observations_ij = [observations_matrix[i, j] for observations_matrix in list_of_observations]
                sample_variance = statistics.variance(observations_ij)
                variance = variance_matrix[i, j]

                sample_variances.append(sample_variance)
                variances.append(variance)

        sample_variances = [np.log(np.log(s)) for s in sample_variances if s >= 1.0]
        variances = [np.log(np.log(s)) for s in variances if s >= 1.0]

        print('Plot')
        plt.subplot(2, 2, plot_num)
        plt.plot([0.0, np.log(np.log(1000.0))], [0.0, np.log(np.log(1000.0))])

        # scatter plot
        plt.scatter(variances, sample_variances, s=2, color='hotpink')
        plt.title(f'alpha: {alpha}. uncertainty: {uncertainty}')
        plt.xlabel('Variance from formula')
        plt.ylabel('Variance over experiments')
plt.show()
