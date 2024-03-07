import random
import numpy as np
import scipy.stats as stats


#####
# Special comments:
# EVEN: means we take care of the case i!=j different from the case i=j
# DEBUG: didn't debug this part yet
#####


# need to calculate only on the upper triangle because the matrices are symmetric
def calculate_chi_squared_value(population_amount_, counts_expected_, observed_probabilities_, correction_,
                                should_use_new_test_, cutoff):
    value = 0
    amount_of_small_expected_ = 0
    for row in range(counts_expected_.shape[0]):
        for col in range(row, counts_expected_.shape[1]):
            expected_ = counts_expected_[row, col]
            observed_ = population_amount_ * observed_probabilities_[row, col]
            variance_ = expected_
            if should_use_new_test_:
                if expected_ * correction_[row, col] < cutoff:
                    amount_of_small_expected_ += 1
                    continue
            else:
                if expected_ < cutoff:
                    amount_of_small_expected_ += 1
                    continue
            if should_use_new_test_:
                # value += correction_[row, col] * (((expected_ - observed_) ** 2) / variance_)
                value += (1 / correction_[row, col]) * (((expected_ - observed_) ** 2) / variance_)
            else:
                value += (((expected_ - observed_) ** 2) / variance_)
    return value, amount_of_small_expected_


def run_experiment(alleles_count, population_amount, alleles_probabilities,
                   observed_probabilities, correction, should_use_new_test, cutoff_value_):
    counts_expected = np.zeros((alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # EVEN
            mult = 1
            if k != j:
                mult = 2
            expected_value = mult * population_amount * alleles_probabilities[k] * alleles_probabilities[j]
            counts_expected[k, j] = expected_value
            counts_expected[j, k] = expected_value

    chi_squared_stat, amount_of_small_expected = calculate_chi_squared_value(population_amount_=population_amount,
                                                                             counts_expected_=counts_expected,
                                                                             observed_probabilities_=observed_probabilities,
                                                                             correction_=correction,
                                                                             should_use_new_test_=should_use_new_test,
                                                                             cutoff=cutoff_value_)
    couples_amount = (alleles_count * (alleles_count + 1)) / 2 - 1
    dof = couples_amount - amount_of_small_expected

    # print(f' alpha for choice: {alpha_val}')
    # print(f' chi square value: {chi_squared_stat}')

    # crit = stats.chi2.ppf(q=0.95, df=dof)
    # print(f'Critical value: {crit}')

    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,
                                 df=dof)
    return p_value, chi_squared_stat, dof