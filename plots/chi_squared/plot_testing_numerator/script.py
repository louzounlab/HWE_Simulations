import random
import numpy as np
import math
import statistics
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import utils_with_certainty


def prepare_probabilities(alleles_count):
    # probabilities {p(i)}
    alleles_probabilities = utils_with_certainty.calculate_alleles_probabilities(alleles_count)

    return alleles_probabilities


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   alleles_probabilities):
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
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

    # now we multiply the upper triangle by 2
    for t in range(alleles_count):
        for m in range(t + 1, alleles_count):
            probabilities[t, m] *= 2

    # print(f'alleles_individuals: {alleles_individuals}')

    # # count the amounts of alleles occurrences in the population
    counts_observed = np.zeros((alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # count amount of occurrences of j,k
            count_tuple = utils_with_certainty.count_row_occurrence_in_2d_array(j, k, alleles_individuals)
            counts_observed[j, k] = count_tuple

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
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            probability = 0.0
            for k in range(population_amount):
                probability += all_probabilities[t, m, k]
            probability /= population_amount

            observed_probabilities[t, m] = probability

    observations = observed_probabilities * population_amount

    # every element is (O(i,j) - P_a(i,j)) ^ 2
    observations_distances_list = []
    # every element is (counts_observed(i,j) - P_a(i,j)) ^ 2
    counts_observed_distances_list = []

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            observed_distance = (observations[i, j] - population_amount * probabilities[i, j]) ** 2
            counts_observed_distance = (counts_observed[i, j] - population_amount * probabilities[i, j]) ** 2

            observations_distances_list.append(observed_distance)
            counts_observed_distances_list.append(counts_observed_distance)

    return observations_distances_list, counts_observed_distances_list


alleles_amount = 20  # 20
population_size = 10000  # 10000
alpha_vals = [1.0, 0.8]
uncertainty_vals = [0.0, 0.2, 0.4]

plt.rcParams["font.family"] = "Arial"
fig, axes = plt.subplots(len(alpha_vals), len(uncertainty_vals), constrained_layout=True)

alleles_probabilities_ = prepare_probabilities(alleles_count=alleles_amount)

plot_num = 0

for alpha in alpha_vals:
    for uncertainty in uncertainty_vals:
        print(f'alpha: {alpha}. uncertainty: {uncertainty}')
        plot_num += 1
        observations_distances_list_, counts_observed_distances_list_ = run_experiment(alleles_count=alleles_amount,
                                                                                       population_amount=population_size,
                                                                                       alpha_val=alpha,
                                                                                       uncertainty_val=uncertainty,
                                                                                       alleles_probabilities=alleles_probabilities_)
        observations_distances_list_ = [np.log(np.log(s)) for s in observations_distances_list_]
        counts_observed_distances_list_ = [np.log(np.log(s)) for s in counts_observed_distances_list_]

        plt.subplot(len(alpha_vals), len(uncertainty_vals), plot_num)
        plt.plot([0.0, np.log(np.log(100000.0))], [0.0, np.log(np.log(100000.0))])
        # plt.plot([0.0, 1000.0], [0.0, 1000.0])
        # scatter plot
        plt.scatter(counts_observed_distances_list_, observations_distances_list_, s=2, color='hotpink')
        plt.title(f'alpha: {alpha}. uncertainty: {uncertainty}')
        plt.xlabel('(Z_ij - N * P_alpha(i,j)) ^ 2')
        plt.ylabel('(O_ij - N * P_alpha(i,j)) ^ 2')
plt.show()
