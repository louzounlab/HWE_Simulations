import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import utils
import warnings
import chi_squared_comparison_real_data
from models.gibbs_sampling import full_algorithm
import single_experiment_with_certainty_test

real_data_path = '../../../data/Snipping_Data/snps_from_martin/donors.ped'


# for a given snipping, check that we have 2 alleles. if yes return 1. otherwise return 0
def is_valid_snipping(snipping_index):
    starting_index = 6
    forbidden_char = '0'
    set_letters = set()

    # here we only map the alleles into indices, count the amount of alleles and population
    with open(real_data_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # print(f'row: {index}')
            lst = line.split(' ')
            alleles = lst[starting_index + 2 * snipping_index], lst[starting_index + 2 * snipping_index + 1]

            if forbidden_char in alleles:
                continue
            set_letters.add(alleles[0])
            set_letters.add(alleles[1])
    if len(set_letters) == 2:
        return 1
    return 0


# get the first amount-sized list of valid snippings, given that there are total amount_in_files snippings
def get_first_valid_snippings(amount, amount_in_file):
    valid_snippings = []
    for snipping_index_ in range(amount_in_file):
        if is_valid_snipping(snipping_index_):
            valid_snippings.append(snipping_index_)
        if len(valid_snippings) == amount:
            break
    return valid_snippings


# given a snipping index. get all the observations of the donors for that snipping.
# return: alleles_count, population_amount, alleles_probabilities, observed_probabilities, observations, probabilities
def get_data(snipping_index):
    # from this index all the snipping appear
    starting_index = 6

    # dictionary: {G -> 0, A -> 1, ...}
    alleles_to_indices = {}

    population_amount = 0

    forbidden_char = '0'

    # here we only map the alleles into indices, count the amount of alleles and population
    with open("../../../data/Snipping_Data/snps_from_martin/donors.ped", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # print(f'row: {index}')
            lst = line.split(' ')
            alleles = lst[starting_index + 2 * snipping_index], lst[starting_index + 2 * snipping_index + 1]

            # if one of the alleles is '0'. its missing value and we need to skip.
            if (alleles[0] == forbidden_char) or (alleles[1] == forbidden_char):
                continue
            else:
                population_amount += 1

            if alleles[0] not in alleles_to_indices:
                alleles_to_indices[alleles[0]] = len(alleles_to_indices)
            if alleles[1] not in alleles_to_indices:
                alleles_to_indices[alleles[1]] = len(alleles_to_indices)

    alleles_count = len(alleles_to_indices)
    observations = np.zeros(shape=(alleles_count, alleles_count))

    # here we get all the observations
    with open("../../../data/Snipping_Data/snps_from_martin/donors.ped", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # print(f'row: {index} / {population_amount}')
            lst = line.split(' ')
            # get the two alleles from the snipping
            alleles = lst[starting_index + 2 * snipping_index], lst[starting_index + 2 * snipping_index + 1]
            # if one of the alleles is '0'. its missing value and we need to skip.
            if (alleles[0] == forbidden_char) or (alleles[1] == forbidden_char):
                continue
            # get the indices of the alleles
            i, j = [alleles_to_indices[allele] for allele in alleles]
            # make sure i <= j
            i, j = min(i, j), max(i, j)

            observations[i, j] += 1

    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            observed_probabilities[i, j] = observations[i, j] / population_amount

    alleles_probabilities = np.zeros(alleles_count)

    for i in range(alleles_count):
        probability = 0.0
        for j in range(alleles_count):
            if i == j:
                probability += observed_probabilities[i, i]
            else:
                probability += observed_probabilities[min(i, j), max(i, j)] * 0.5
        alleles_probabilities[i] = probability

    # print('checking summ is 1:')
    #
    # sum = 0.0
    # for i in range(alleles_count):
    #     sum += alleles_probabilities[i]
    # print(f'alleles probabilities: {sum}')
    #
    # sum = 0.0
    # for i in range(alleles_count):
    #     for j in range(i, alleles_count):
    #         sum += probabilities[i, j]
    # print(f'probabilities: {sum}')

    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            mult = 2.0
            if i == j:
                mult = 1.0
            probabilities[i, j] = mult * alleles_probabilities[i] * alleles_probabilities[j]

    return alleles_count, population_amount, alleles_probabilities, observed_probabilities, observations, probabilities


if __name__ == '__main__':
    warnings.filterwarnings("error")
    # for every snipping index: perform gibbs sampling and chi squared.
    # make a scatter plot: (index, result).
    gibbs_sampling_results = []
    chi_squared_results = []

    # got from the test file: get length l of a row and calc: (l - 6) / 2
    snipping_amount_from_file = 13258
    snipping_amount = 100
    #
    print('Getting valid snippings')
    snipping_values = get_first_valid_snippings(snipping_amount, snipping_amount_from_file)
    print('Got valid snippings')

    for i, snipping_index_ in enumerate(snipping_values):
        print(f'current index: {i} / {snipping_amount - 1}, current snipping: {snipping_index_}')
        print('getting observations')
        alleles_count_, population_amount_, \
            alleles_probabilities_, \
            observed_probabilities_, \
            observations_, \
            probabilities_ = get_data(snipping_index_)
        observed_cdf_ = single_experiment_with_certainty_test.calc_observed_cdf(alleles_count=alleles_count_,
                                                                                observed=observations_)

        print('getting chi-squared p_value')
        # get the chi-squared p_value
        chi_squared_p_value, _, _ = chi_squared_comparison_real_data.run_experiment(alleles_count=alleles_count_,
                                                                                    population_amount=population_amount_,
                                                                                    uncertainty_val=0.0,
                                                                                    alleles_probabilities=alleles_probabilities_,
                                                                                    observed_probabilities=observed_probabilities_)
        print('getting gibbs sampling p_value')
        # get the gibbs sampling p_value
        gibbs_sampling_p_value = full_algorithm(alleles_count=alleles_count_,
                                                population_amount=population_amount_,
                                                alleles_probabilities=alleles_probabilities_,
                                                probabilities=probabilities_,
                                                observations=observations_,
                                                observed_cdf=observed_cdf_)

        chi_squared_results.append(chi_squared_p_value)
        gibbs_sampling_results.append(gibbs_sampling_p_value)

    columns = ['snipping_index', 'gibbs_sampling_p_value', 'chi_squared_p_value']
    index = list(range(snipping_amount))
    df = pd.DataFrame(index=index, columns=columns)
    df['snipping_index'] = snipping_values
    df['gibbs_sampling_p_value'] = gibbs_sampling_results
    df['chi_squared_p_value'] = chi_squared_results
    df.to_csv('results.csv')
