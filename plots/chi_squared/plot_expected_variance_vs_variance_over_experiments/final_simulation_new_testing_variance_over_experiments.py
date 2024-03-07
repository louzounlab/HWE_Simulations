import random
import numpy as np
import math
import statistics
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from models_simulated import chi_squared


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

    # print(f'alleles: {alleles_probabilities}')
    # print(f' marginal: {marginal_probabilities}')

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

    # matrix p_k(i,j) = A[i,j,k]
    all_probabilities = np.zeros(shape=(alleles_count, alleles_count, population_amount))

    # print(f'alleles_individuals after uncertainty: {alleles_individuals}')

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
    return observations


if __name__ == '__main__':
    num_of_experiments = 100  # 500 amount of experiments for each alpha ranging from zero to one
    alleles_amount = 20  # 20
    population_size = 10000  # 10000
    alpha_vals = [1.0, 0.8]
    uncertainty_vals = [0.0, 0.2, 0.4]
    # alpha_vals = [1.0]
    # uncertainty_vals = [0.0]

    sns.set_style('white')
    plt.rcParams["font.family"] = "Arial"

    fig, axes = plt.subplots(len(alpha_vals), len(uncertainty_vals), constrained_layout=True)

    print('Calculating probabilities')
    # get {p(i)}, {p(i|j)}
    alleles_probabilities_ = chi_squared.prepare_probabilities(alleles_count=alleles_amount)
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
                observations_ = run_experiment(alleles_count=alleles_amount,
                                               population_amount=population_size,
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
            min_val = min(variances_ + sample_variances_)
            max_val = max(variances_ + sample_variances_)
            plt.xlim(min_val, max_val)
            plt.ylim(min_val, max_val)
            ax = plt.gca()
    plt.savefig('plot.png', pad_inches=0.2, bbox_inches="tight")
    plt.show()
