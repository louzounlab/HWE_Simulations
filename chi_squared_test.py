import numpy as np
import pandas as pd
import scipy.stats as stats
import utils


def calculate_chi_squared_value(alleles_amount_, population_amount_, alleles_probabilities_, observed_,
                                counts_total_variance_, should_use_new_test_):
    value = 0.0
    for i in range(alleles_amount_):
        for j in range(i, alleles_amount_):
            mult = 1
            if i != j:
                mult = 2
            expected_val = population_amount_ * alleles_probabilities_[i] * alleles_probabilities_[j] * mult
            observed_val = observed_[i, j]
            variance_val = counts_total_variance_[i, j]
            if should_use_new_test_:
                denominator = variance_val
            else:
                denominator = expected_val
            value += (((expected_val - observed_val) ** 2) / denominator)
            # print(f'expected value: {expected_val}, observed_val : {observed_val}, denominator: {denominator}')
    return value


def run_experiment(alleles_count, population_amount, uncertainty_val, alleles_probabilities,
                   observed,
                   should_use_new_test_):
    # i_j -> {i,j,V(i,j)}
    total_variance_matrix = np.zeros(shape=(alleles_count, alleles_count))
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            # EVEN
            mult = 1
            if i != j:
                mult = 2
            expected = population_amount * alleles_probabilities[i] * alleles_probabilities[j] * mult
            total_variance_matrix[i, j] = expected + observed[i, j] * uncertainty_val

    chi_squared_stat = calculate_chi_squared_value(alleles_amount_=alleles_count,
                                                   population_amount_=population_amount,
                                                   alleles_probabilities_=alleles_probabilities,
                                                   observed_=observed,
                                                   counts_total_variance_=total_variance_matrix,
                                                   should_use_new_test_=should_use_new_test_)
    print(f' statistic: {chi_squared_stat}')
    dof = (alleles_count * (alleles_count + 1)) / 2 - 1
    print(f'dof: {dof}')

    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,
                                 df=dof)

    print(f'p_value: {p_value}')


alleles_amount = 10
alleles_probabilities_ = np.zeros(10) + 0.1
population_amount_ = 1000
uncertainty = 0.0
observed_ = np.zeros(shape=(10, 10))
for i in range(10):
    for j in range(i, 10):
        if i == j:
            observed_[i, j] = population_amount_ * (alleles_probabilities_[i] ** 2)
        else:
            observed_[i, j] = 2 * population_amount_ * alleles_probabilities_[i] * alleles_probabilities_[j]
observed_ += 5
should_use_new_test = 0

run_experiment(alleles_count=alleles_amount, population_amount=population_amount_,
               uncertainty_val=uncertainty, alleles_probabilities=alleles_probabilities_,
               observed=observed_, should_use_new_test_=should_use_new_test)
