import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import utils
from models import chi_squared

real_data_path = '../../../data/hla_data_5_populations'


def get_data_for_chi_squared():
    levels_set = utils.LEVELS_LIST
    races_set = utils.RACES_5_LIST

    columns = list(races_set)
    index = list(levels_set)

    df_old_chi_squared_results = pd.DataFrame(columns=columns, index=index)
    df_new_chi_squared_results = pd.DataFrame(columns=columns, index=index)
    df_sampling_chi_squared_results = pd.DataFrame(columns=columns, index=index)
    df_old_chi_squared = pd.DataFrame(columns=columns, index=index)
    df_new_chi_squared = pd.DataFrame(columns=columns, index=index)
    df_sampling_chi_squared = pd.DataFrame(columns=columns, index=index)
    # df_ambiguities = pd.DataFrame(columns=columns, index=index)
    df_dof = pd.DataFrame(columns=columns, index=index)

    for i, level in enumerate(index):
        for j, race in enumerate(columns):
            print(f'level: {level}, race: {race}')
            old_p_value_, new_p_value_, sampling_p_value_,\
                old_chi_squared_, new_chi_squared_, sampling_chi_squared_, \
                dof_ = perform_tests_save_data(level, race)

            df_old_chi_squared_results.iloc[i, j] = old_p_value_
            df_new_chi_squared_results.iloc[i, j] = new_p_value_
            df_sampling_chi_squared_results.iloc[i, j] = sampling_p_value_

            df_old_chi_squared.iloc[i, j] = old_chi_squared_
            df_new_chi_squared.iloc[i, j] = new_chi_squared_
            df_sampling_chi_squared.iloc[i, j] = sampling_chi_squared_

            # df_ambiguities.iloc[i, j] = ambiguity
            df_dof.iloc[i, j] = dof_

    df_old_chi_squared_results.to_csv(f'{real_data_path}/results/old_chi_squared_p_values')
    df_new_chi_squared_results.to_csv(f'{real_data_path}/results/new_chi_squared_p_values')
    df_sampling_chi_squared_results.to_csv(f'{real_data_path}/results/sampling_chi_squared_p_values')

    df_old_chi_squared.to_csv(f'{real_data_path}/results/old_chi_squared_statistic')
    df_new_chi_squared.to_csv(f'{real_data_path}/results/new_chi_squared_statistic')
    df_sampling_chi_squared.to_csv(f'{real_data_path}/results/sampling_chi_squared_statistic')
    # df_ambiguities.to_csv(f'{real_data_path}/results/ambiguities')
    df_dof.to_csv(f'{real_data_path}/results/dof')


def perform_tests_save_data(level, race):
    # df_id_allele1_allele2_prob = pd.read_csv(
    #     f'{real_data_path}/levels/{level}/races/{race}/id_allele1_allele2_probability',
    #     dtype={'id': str, 'allele_1': int, 'allele_2': int, 'probability': float})
    id_to_index = {}
    allele_to_index = {}

    probabilities_path = f'{real_data_path}/levels/{level}/races/{race}/id_allele1_allele2_probability'

    # first read all the rows and get indices of ids and alleles and amounts
    with open(probabilities_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # header
            if index == 0:
                continue
            if index % 10000 == 0:
                print(f'row: {index}, level: {level}, race: {race}, preprocessing')
            lst = line.strip('\n').split(',')

            id = lst[0]
            allele_1 = lst[1]
            allele_2 = lst[2]
            probability = float(lst[3])

            if id not in id_to_index:
                id_to_index[id] = len(id_to_index)

            if allele_1 not in allele_to_index:
                allele_to_index[allele_1] = len(allele_to_index)
            if allele_2 not in allele_to_index:
                allele_to_index[allele_2] = len(allele_to_index)

    # first read all the rows and get indices of ids and alleles and amounts
    # for index, row in df_id_allele1_allele2_prob.iterrows():
    #     if index % 10000 == 0:
    #         print(f'row: {index}, level: {level}, race: {race}, preprocessing')
    #     id = row['id']
    #     allele_1 = row['allele_1']
    #     allele_2 = row['allele_2']
    #     probability = row['probability']
    #
    #     if id not in id_to_index:
    #         id_to_index[id] = len(id_to_index)
    #
    #     if allele_1 not in allele_to_index:
    #         allele_to_index[allele_1] = len(allele_to_index)
    #     if allele_2 not in allele_to_index:
    #         allele_to_index[allele_2] = len(allele_to_index)

    alleles_count = len(allele_to_index)
    population_amount = len(id_to_index)

    # {p(i)}
    alleles_probabilities = np.zeros(alleles_count)

    # {p(i, j)}
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))

    # correction matrix
    correction = np.zeros(shape=(alleles_count, alleles_count))

    # calculate {p_k(i,j)}
    with open(probabilities_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index == 0:
                continue

            if index % 10000 == 0:
                print(f'row: {index}, level: {level}, race: {race}, After preprocessing')
            lst = line.strip('\n').split(',')

            # id = lst[0]
            allele_1 = lst[1]
            allele_2 = lst[2]
            probability = float(lst[3])

            # id_index = id_to_index[id]

            allele_1_index = allele_to_index[allele_1]
            allele_2_index = allele_to_index[allele_2]

            allele_1_index, allele_2_index = min(allele_1_index, allele_2_index), max(allele_1_index, allele_2_index)

            alleles_probabilities[allele_1_index] += 0.5 * probability
            alleles_probabilities[allele_2_index] += 0.5 * probability

            observed_probabilities[allele_1_index, allele_2_index] += probability

            correction[allele_1_index, allele_2_index] += (probability ** 2)

    alleles_probabilities /= population_amount

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            observed_probabilities[i, j] /= population_amount
            if observed_probabilities[i, j] == 0:
                correction[i, j] = 1.0
            else:
                correction[i, j] /= (population_amount * observed_probabilities[i, j])

    print('old test')
    old_p_value, old_chi_squared, dof = chi_squared.run_experiment(alleles_count=alleles_count,
                                                                   population_amount=population_amount,
                                                                   alleles_probabilities=alleles_probabilities,
                                                                   observed_probabilities=observed_probabilities,
                                                                   correction=correction,
                                                                   should_use_new_test=0,
                                                                   cutoff_value_=2.0)
    print(f'p_value: {old_p_value}')

    print('new test')
    new_p_value, new_chi_squared, _ = chi_squared.run_experiment(alleles_count=alleles_count,
                                                                 population_amount=population_amount,
                                                                 alleles_probabilities=alleles_probabilities,
                                                                 observed_probabilities=observed_probabilities,
                                                                 correction=correction,
                                                                 should_use_new_test=1,
                                                                 cutoff_value_=2.0)
    print(f'p_value: {new_p_value}')

    # here we also perform Chi Squared sampling (first sample from {p_k(i,j)}
    # the ids are sorted, for each donor k we sample alleles (i, j) from his set of probabilities
    observed_probabilities_sampling = np.zeros(shape=(alleles_count, alleles_count))
    # the last id that we scanned
    last_id = ''
    # the probabilities of the last id. list of floats
    last_probabilities = []
    # the possible alleles of the last id, according to his probabilities. list where each element is a list [i, j]
    last_possible_alleles = []

    with open(probabilities_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index == 0:
                continue
            lst = line.strip('\n').split(',')
            # '47919664'
            id_row = lst[0]

            # if we finished scanning the last id
            if last_id and (id_row != last_id):
                # sample alleles for the last id
                index = random.choices(population=range(len(last_probabilities)), weights=last_probabilities,
                                       k=1)[0]
                # get the sampled alleles
                i, j = last_possible_alleles[index]
                # update observations
                observed_probabilities_sampling[i, j] += 1
                # reset last person
                last_probabilities = []
                last_possible_alleles = []

            allele_1 = lst[1]
            allele_2 = lst[2]
            probability = float(lst[3])

            # id_index = id_to_index[id]

            allele_1_index = allele_to_index[allele_1]
            allele_2_index = allele_to_index[allele_2]

            allele_1_index, allele_2_index = min(allele_1_index, allele_2_index), max(allele_1_index, allele_2_index)

            last_id = id_row
            last_probabilities.append(probability)
            last_possible_alleles.append([allele_1_index, allele_2_index])

    # we still have the last id
    index = random.choices(population=range(len(last_probabilities)), weights=last_probabilities,
                           k=1)[0]
    # get the sampled alleles
    i, j = last_possible_alleles[index]
    # update observations
    observed_probabilities_sampling[i, j] += 1

    # normalize observations into probabilities
    observed_probabilities_sampling /= np.sum(observed_probabilities_sampling)

    alleles_probabilities_sampling = np.zeros(alleles_count)

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            alleles_probabilities_sampling[i] += 0.5 * observed_probabilities_sampling[i, j]
            alleles_probabilities_sampling[j] += 0.5 * observed_probabilities_sampling[i, j]

    print('sampling test')
    sampling_p_value, sampling_chi_squared, _ = chi_squared.run_experiment(alleles_count=alleles_count,
                                                                           population_amount=population_amount,
                                                                           alleles_probabilities=alleles_probabilities_sampling,
                                                                           observed_probabilities=observed_probabilities_sampling,
                                                                           correction=correction,
                                                                           should_use_new_test=0,
                                                                           cutoff_value_=2.0)
    print(f'sampling p value: {sampling_p_value}')

    return old_p_value, new_p_value, sampling_p_value, old_chi_squared, new_chi_squared, sampling_chi_squared, dof


if __name__ == '__main__':
    get_data_for_chi_squared()
