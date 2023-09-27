import random
import numpy as np
import math
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import utils_with_uncertainty


# This simulation uses matrices instead of dictionaries
def run_experiment(alleles_count, population_amount, uncertainty_val, alleles_probabilities,
                   observed,
                   should_use_new_test_):
    counts_total_variance = np.zeros((alleles_count, alleles_count))
    for k in range(alleles_count):
        for j in range(k, alleles_count):
            total_variance = utils_with_uncertainty.calculate_total_variance_alleles_i_j(alleles_count,
                                                                                         population_amount,
                                                                                         alleles_probabilities,
                                                                                         uncertainty_val,
                                                                                         observed,
                                                                                         k, j)

            counts_total_variance[k, j] = total_variance
            counts_total_variance[j, k] = total_variance

    chi_squared_stat = utils_with_uncertainty.calculate_chi_squared_value(alleles_amount=alleles_count,
                                                                          population_amount=population_amount,
                                                                          alleles_probabilities=alleles_probabilities,
                                                                          counts_observed_=observed,
                                                                          counts_total_variance_=counts_total_variance,
                                                                          should_use_new_test_=should_use_new_test_)
    dof = (alleles_count * (alleles_count + 1)) / 2 - 1

    # print(f' alpha for choice: {alpha_val}')
    # print(f' chi square value: {chi_squared_stat}')

    # crit = stats.chi2.ppf(q=0.95, df=dof)
    # print(f'Critical value: {crit}')

    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,
                                 df=dof)
    # if p_value < 0.05:
    #     return 1
    # return 0

    return p_value
