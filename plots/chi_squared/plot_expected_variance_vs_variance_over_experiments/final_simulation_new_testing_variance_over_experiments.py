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
    sum_ = 0.0
    for row in range(x.shape[0]):
        sum_ += np.exp(x[row])
    return np.exp(x) / sum_


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

# given a matrix of probabilities. take the cdf of the upper triangle part
def get_cdf_list_from_matrix(alleles_count, probabilities):
    couples = alleles_count * (alleles_count + 1) // 2
    cdf = np.zeros(couples)
    index = 0
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            cdf[index] = probabilities[i, j]
            if index > 0:
                cdf[index] += cdf[index - 1]
            index += 1
    return cdf

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

    return alleles_probabilities


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   alleles_probabilities):

    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
                probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
    # matrix where every row is the alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, 2), dtype=np.int32)

    # 1) assignment of alleles with certainty
    for k in range(population_amount):
        probabilities_list = probabilities.flatten()

        # notice the alleles i != j have probability 2 * p(i,j)
        index = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=1)[0]
        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row , col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

    # now we multiply the upper triangle by 2
    for t in range(alleles_count):
        for m in range(t + 1, alleles_count):
            probabilities[t, m] *= 2

    # print(f'alleles_individuals: {alleles_individuals}')

    # # count the amounts of alleles occurrences in the population
    # counts_observed = np.zeros((alleles_count, alleles_count))
    # for j in range(alleles_count):
    #     for k in range(j, alleles_count):
    #         # count amount of occurrences of j,k
    #         count_tuple = count_row_occurrence_in_2d_array(j, k, alleles_individuals)
    #
    #         counts_observed[j, k] = count_tuple

    # matrix p_k(i,j) = A[i,j,k]
    all_probabilities = np.zeros(shape=(alleles_count, alleles_count, population_amount))

    # adding uncertainty to our model

    for k in range(population_amount):
        # person k has alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)

        # choice whether this person will have uncertain alleles
        choice = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=1)[0]

        # this person has certain alleles
        if choice == 1:
            all_probabilities[j, l, k] = 1.0
        # this person has uncertain alleles
        if choice == 0:
            for t in range(alleles_count):
                for m in range(t, alleles_count):
                    all_probabilities[t, m, k] = probabilities[t, m]

    # now we have the probabilities {p_k(i,j)}
    # lets get the final p(i,j) = 1 / N * sum_k p_k(i,j)

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            probability = 0.0
            for k in range(population_amount):
                probability += all_probabilities[t, m, k]
            probability /= population_amount

            probabilities[t, m] = probability

    observations = probabilities * population_amount
    return observations



num_of_experiments = 100  # 500 amount of experiments for each alpha ranging from zero to one
alleles_amount = 20  # 20
population_size = 10000  # 10000
alpha_vals = [1.0, 0.8]
uncertainty_vals = [0.0, 0.2, 0.4]
# alpha_vals = [1.0]
# uncertainty_vals = [0.0]


plt.rcParams["font.family"] = "Arial"
fig, axes = plt.subplots(len(alpha_vals), len(uncertainty_vals), constrained_layout=True)

print('Calculating probabilities')
# get {p(i)}, {p(i|j)}
alleles_probabilities_ = prepare_probabilities(alleles_count=alleles_amount)
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
                variance = mult * population_size * alleles_probabilities_[i] * alleles_probabilities_[j] * \
                           (1 + uncertainty)
                variance_matrix[i, j] = variance

        # define a list where every element is a matrix {O(i,j)}
        list_of_observations = []

        print('Running experiments')

        for experiment in range(num_of_experiments):
            print(
                f'Experiment num: {experiment} / {num_of_experiments}: plot: {plot_num} / {len(alpha_vals) * len(uncertainty_vals)}')
            observations_ = run_experiment(alleles_count=alleles_amount, population_amount=population_size,
                                           alpha_val=alpha, uncertainty_val=uncertainty,
                                           alleles_probabilities=alleles_probabilities_)
            list_of_observations.append(observations_)

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

        # sample_variances_ = [np.log(np.log(s)) for s in sample_variances if s > 0.0]
        # variances_ = [np.log(np.log(variances[i])) for i in range(len(variances)) if sample_variances[i] > 0.0]
        sample_variances_ = [np.log(np.log(s)) for s in sample_variances]
        variances_ = [np.log(np.log(variances[i])) for i in range(len(variances))]

        print('Plot')
        plt.subplot(len(alpha_vals), len(uncertainty_vals), plot_num)
        plt.plot([0.0, np.log(np.log(1000.0))], [0.0, np.log(np.log(1000.0))])

        # scatter plot
        plt.scatter(variances_, sample_variances_, s=2, color='hotpink')
        plt.title(f'alpha: {alpha}. uncertainty: {uncertainty}')
        plt.xlabel('Variance from formula')
        plt.ylabel('Variance over experiments')
plt.show()
