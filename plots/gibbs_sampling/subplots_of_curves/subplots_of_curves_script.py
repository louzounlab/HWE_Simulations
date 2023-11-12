import numpy as np
import matplotlib.pyplot as plt
from models_simulated import gibbs_sampling


def subplots_of_curves(alpha_val):
    alleles_count_ = 1000
    population_amount_ = 1000000

    observed_, _ \
        = gibbs_sampling.prepare_experiment_data(alleles_count=alleles_count_,
                                                 population_amount=population_amount_,
                                                 alpha_val=alpha_val)
    population_amount_ = np.sum(observed_)  # check this is integer

    # in case the matrix is not upper triangular
    for i in range(alleles_count_):
        for j in range(i + 1, alleles_count_):
            observed_[i, j] += observed_[j, i]

    # observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))
    # for i in range(alleles_count):
    #     for j in range(i, alleles_count):
    #         observed_probabilities[i, j] = observations[i, j] / population_amount_calculated

    alleles_probabilities = np.zeros(alleles_count_)

    for i in range(alleles_count_):
        probability = 0.0
        for j in range(alleles_count_):
            if i == j:
                probability += observed_[i, i]
            else:
                probability += observed_[min(i, j), max(i, j)] * 0.5
        alleles_probabilities[i] = probability / population_amount_

    probabilities = np.zeros(shape=(alleles_count_, alleles_count_))
    for i in range(alleles_count_):
        for j in range(i, alleles_count_):
            mult = 2.0
            if i == j:
                mult = 1.0
            probabilities[i, j] = mult * alleles_probabilities[i] * alleles_probabilities[j]

    observed_cdf_ = gibbs_sampling.calc_observed_cdf(alleles_count_, observed_)

    observed_copy = np.copy(observed_)
    observed_cdf_copy = np.copy(observed_cdf_)

    # plt.title(title, weight='bold', fontsize=20)
    plt.subplot(7, 7, 1)
    for i in range(1, 50):
        print(f'alpha: {alpha_val}, plot: {i} / 49')

        gibbs_sampling.perform_experiment(alleles_count_, population_amount_, alleles_probabilities,
                                          probabilities,
                                          observed_copy, observed_cdf_copy, i)
        np.copyto(observed_copy, observed_)
        np.copyto(observed_cdf_copy, observed_cdf_)
    alpha_str = str(int(alpha_val * 1000))
    plt.savefig(f'subplots_of_curves_alpha_{alpha_str}', pad_inches=0.2, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    alpha_vals = [1.0, 0.99]

    for i, alpha_val_ in enumerate(alpha_vals):
        subplots_of_curves(alpha_val=alpha_val_)
