import numpy as np
import pandas as pd
import scipy.stats as stats
import utils

real_data_path = 'data/real_data_test'


def calculate_total_variance_alleles_i_j(alleles_amount_,
                                         population_amount,
                                         alleles_probabilities,
                                         uncertainty_,
                                         observed_,
                                         i_, j_):
    # EVEN
    mult = 1
    if i_ != j_:
        mult = 2
    expected = population_amount * alleles_probabilities[i_] * alleles_probabilities[j_] * mult

    # i_j_observed_probability = observed_[i_, j_] / sum_observed_
    i_j_observed_probability = observed_[i_, j_]

    return expected + i_j_observed_probability * uncertainty_


def calculate_chi_squared_value(alleles_amount, population_amount, alleles_probabilities, observed_probabilities_,
                                counts_total_variance_, should_use_new_test_):
    value = 0.0
    small_expected_observed_counter = 0
    small_expected_high_observed_counter = 0
    for i in range(alleles_amount):
        for j in range(i, alleles_amount):
            mult = 1
            if i != j:
                mult = 2
            expected_val = population_amount * alleles_probabilities[i] * alleles_probabilities[j] * mult
            observed_probability = observed_probabilities_[i, j]
            observed_val = population_amount * observed_probability
            # print(f'expected: {expected_val}')
            # print(f'observed: {observed_val}')
            if (expected_val < 2) and (observed_val < 2):
                small_expected_observed_counter += 1
                continue
            if (expected_val < 2) and (observed_val > 2):
                small_expected_high_observed_counter += 1
            variance_val = counts_total_variance_[i, j]
            if should_use_new_test_:
                denominator = variance_val
            else:
                denominator = expected_val
            value += (((expected_val - observed_val) ** 2) / denominator)
            # print(f'expected value: {expected_val}, observed_val : {observed_val}, denominator: {denominator}')
    print(f'amount of small expected and high observed is: {small_expected_high_observed_counter}')
    return value, small_expected_observed_counter


# alleles_probabilities, observed_ are dictionaries
def run_experiment(alleles_count, population_amount, uncertainty_val, alleles_probabilities,
                   observed_probabilities,
                   should_use_new_test_):
    # i_j -> {i,j,V(i,j)}
    total_variance_matrix = np.zeros(shape=(alleles_count, alleles_count))
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            # EVEN
            mult = 1
            if i != j:
                mult = 2
            expected_probability = alleles_probabilities[i] * alleles_probabilities[j] * mult
            # this is equivalent to: N*p(i)p(j) + N*p(i)p(j)*u where u is the average uncertainty
            total_variance_matrix[i, j] = population_amount * expected_probability + expected_probability * uncertainty_val

    chi_squared_stat, small_expected_small_observed_counter = calculate_chi_squared_value(alleles_amount=alleles_count,
                                                                                          population_amount=population_amount,
                                                                                          alleles_probabilities=alleles_probabilities,
                                                                                          observed_probabilities_=observed_probabilities,
                                                                                          counts_total_variance_=total_variance_matrix,
                                                                                          should_use_new_test_=should_use_new_test_)
    print(f'amount of small expected and observed: {small_expected_small_observed_counter}')
    print(f' statistic: {chi_squared_stat}')
    dof = (alleles_count * (alleles_count + 1)) / 2 - 1 - small_expected_small_observed_counter
    print(f'dof: {dof}')
    # print(f' alpha for choice: {alpha_val}')
    # print(f' chi square value: {chi_squared_stat}')

    # crit = stats.chi2.ppf(q=0.95, df=dof)
    # print(f'Critical value: {crit}')

    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,
                                 df=dof)
    # if p_value < 0.05:
    #     return 1
    # return 0

    return p_value, chi_squared_stat, dof


def get_data_for_chi_squared():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    columns = list(races_set)
    index = list(levels_set)

    df_old_chi_squared_results = pd.DataFrame(columns=columns, index=index)
    df_new_chi_squared_results = pd.DataFrame(columns=columns, index=index)
    df_old_chi_squared = pd.DataFrame(columns=columns, index=index)
    df_new_chi_squared = pd.DataFrame(columns=columns, index=index)
    df_ambiguities = pd.DataFrame(columns=columns, index=index)
    df_dof = pd.DataFrame(columns=columns, index=index)

    for i, level in enumerate(index):
        for j, race in enumerate(columns):
            print(f'level: {level}, race: {race}')
            old_p_value_, new_p_value_, old_chi_squared_, new_chi_squared_, ambiguity, dof_ = perform_tests_save_data(
                level,
                race)
            df_old_chi_squared_results.iloc[i, j] = old_p_value_
            df_new_chi_squared_results.iloc[i, j] = new_p_value_
            df_old_chi_squared.iloc[i, j] = old_chi_squared_
            df_new_chi_squared.iloc[i, j] = new_chi_squared_
            df_ambiguities.iloc[i, j] = ambiguity
            df_dof.iloc[i ,j] = dof_

    df_old_chi_squared_results.to_csv(f'{real_data_path}/results/old_chi_squared_p_values')
    df_new_chi_squared_results.to_csv(f'{real_data_path}/results/new_chi_squared_p_values')
    df_old_chi_squared.to_csv(f'{real_data_path}/results/old_chi_squared_statistic')
    df_new_chi_squared.to_csv(f'{real_data_path}/results/new_chi_squared_statistic')
    df_ambiguities.to_csv(f'{real_data_path}/results/ambiguities')
    df_dof.to_csv(f'{real_data_path}/results/dof')


def perform_tests_save_data(level, race):
    # get p(i) and alleles amount

    # df of i, p(i)
    i_probability_df = pd.read_csv(f'{real_data_path}/levels/{level}/races/{race}/i_probability',
                                   dtype={'i': int, 'probability': float})
    # i -> p(i)
    alleles_probabilities_dict = {}

    for index, row in i_probability_df.iterrows():
        allele = int(row['i'])
        probability = row['probability']

        if allele not in alleles_probabilities_dict:
            alleles_probabilities_dict[allele] = probability

    # amount of alleles
    alleles_amount = len(alleles_probabilities_dict)

    alleles_probabilities = np.zeros(alleles_amount)
    for i in range(alleles_amount):
        if i in alleles_probabilities_dict:
            alleles_probabilities[i] = alleles_probabilities_dict[i]
        else:
            alleles_probabilities[i] = 0

    # get p(i,j)

    i_j_probability_df = pd.read_csv(f'{real_data_path}/levels/{level}/races/{race}/i_j_probability',
                                     dtype={'i': int, 'j': int, 'probability': float})
    # i_j -> {i,j,p(i,j)} dict from couples of alleles to dictionaries
    i_j_to_i_j_probability_dict = {}

    observed_probabilities = np.zeros(shape=(alleles_amount, alleles_amount))

    for index, row in i_j_probability_df.iterrows():
        i = int(row['i'])
        j = int(row['j'])
        probability = row['probability']

        # only for safety
        i, j = min(i, j), max(i, j)
        i_j = str(i) + '_' + str(j)

        if i_j not in i_j_to_i_j_probability_dict:
            i_j_to_i_j_probability_dict[i_j] = {'i': i, 'j': j, 'probability': probability}
            observed_probabilities[i, j] = probability

    # get ambiguity and population size

    ambiguity_df = pd.read_csv(f'{real_data_path}/levels/{level}/races/{race}/uncertainty',
                               dtype={'id': str, 'uncertainty': float})
    uncertainty = 0.0

    for index, row in ambiguity_df.iterrows():
        uncertainty_row = row['uncertainty']
        uncertainty += uncertainty_row

    population_amount = ambiguity_df.shape[0]


    # expected_matrix = np.zeros(shape=(alleles_amount, alleles_amount))
    # for i in range(alleles_amount):
    #     for j in range(i, alleles_amount):
    #         mult = 1
    #         if i != j:
    #             mult = 2
    #         expected_matrix[i, j] = population_amount * alleles_probabilities[i] * alleles_probabilities[j] * mult
    # df_expected = pd.DataFrame(expected_matrix)
    # df_observed = pd.DataFrame(observed)
    # df_expected.to_csv("expected_matrix")
    # df_observed.to_csv("observed_matrix")

    print('old test')
    old_p_value, old_chi_squared, dof = run_experiment(alleles_count=alleles_amount,
                                                       population_amount=population_amount,
                                                       uncertainty_val=uncertainty,
                                                       alleles_probabilities=alleles_probabilities,
                                                       observed_probabilities=observed_probabilities,
                                                       should_use_new_test_=0)
    print(f'p_value: {old_p_value}')

    print('new test')
    new_p_value, new_chi_squared, _ = run_experiment(alleles_count=alleles_amount,
                                                     population_amount=population_amount,
                                                     uncertainty_val=uncertainty,
                                                     alleles_probabilities=alleles_probabilities,
                                                     observed_probabilities=observed_probabilities,
                                                     should_use_new_test_=1)
    print(f'p_value: {new_p_value}')

    return old_p_value, new_p_value, old_chi_squared, new_chi_squared, uncertainty, dof


get_data_for_chi_squared()
