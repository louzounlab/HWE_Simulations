import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import constants
from models import chi_squared

real_data_path = '../../../data/haplotypes_data'

# THIS FILE USES THE EXTRACTED DATA, PERFORM THE STATISTICAL TESTS AND SAVE IN THE RESULTS DIRECTORY


def get_data_for_chi_squared(file_name, is_using_haplotypes):
    levels_list = constants.LEVELS_LIST
    races_list = constants.FILE_NAME_TO_POPULATIONS[file_name]

    columns = list(races_list)
    if not is_using_haplotypes:
        index = list(levels_list)
    else:
        index = ['USING_HAPLOTYPES']

    df_old_chi_squared_results = pd.DataFrame(columns=columns, index=index)
    df_new_chi_squared_results = pd.DataFrame(columns=columns, index=index)

    df_old_chi_squared = pd.DataFrame(columns=columns, index=index)
    df_new_chi_squared = pd.DataFrame(columns=columns, index=index)
    # df_ambiguities = pd.DataFrame(columns=columns, index=index)
    df_dof = pd.DataFrame(columns=columns, index=index)

    if not is_using_haplotypes:
        for i, level in enumerate(index):
            for j, race in enumerate(columns):
                print(f'level: {level}, race: {race}')
                old_p_value_, new_p_value_,\
                    old_chi_squared_, new_chi_squared_, \
                    dof_ = perform_tests_save_data(file_name=file_name,
                                                   is_using_haplotypes=is_using_haplotypes,
                                                   level_haplotype=level,
                                                   race=race)

                df_old_chi_squared_results.iloc[i, j] = old_p_value_
                df_new_chi_squared_results.iloc[i, j] = new_p_value_

                df_old_chi_squared.iloc[i, j] = old_chi_squared_
                df_new_chi_squared.iloc[i, j] = new_chi_squared_

                # df_ambiguities.iloc[i, j] = ambiguity
                df_dof.iloc[i, j] = dof_
    else:
        for j, race in enumerate(columns):
            print(f'using haplotypes, race: {race}')
            old_p_value_, new_p_value_, \
                old_chi_squared_, new_chi_squared_, \
                dof_ = perform_tests_save_data(file_name=file_name,
                                               is_using_haplotypes=is_using_haplotypes,
                                               level_haplotype=None,
                                               race=race)

            df_old_chi_squared_results.iloc[0, j] = old_p_value_
            df_new_chi_squared_results.iloc[0, j] = new_p_value_

            df_old_chi_squared.iloc[0, j] = old_chi_squared_
            df_new_chi_squared.iloc[0, j] = new_chi_squared_

            # df_ambiguities.iloc[i, j] = ambiguity
            df_dof.iloc[0, j] = dof_

    if not is_using_haplotypes:
        df_old_chi_squared_results.to_csv(f'{real_data_path}/{file_name}/results/levels/old_chi_squared_p_values.csv')
        df_new_chi_squared_results.to_csv(f'{real_data_path}/{file_name}/results/levels/new_chi_squared_p_values.csv')

        df_old_chi_squared.to_csv(f'{real_data_path}/{file_name}/results/levels/old_chi_squared_statistic.csv')
        df_new_chi_squared.to_csv(f'{real_data_path}/{file_name}/results/levels/new_chi_squared_statistic.csv')
        # df_ambiguities.to_csv(f'{real_data_path}/results/ambiguities')
        df_dof.to_csv(f'{real_data_path}/{file_name}/results/levels/dof.csv')
    else:
        df_old_chi_squared_results.to_csv(f'{real_data_path}/{file_name}/results/haplotypes/old_chi_squared_p_values.csv')
        df_new_chi_squared_results.to_csv(f'{real_data_path}/{file_name}/results/haplotypes/new_chi_squared_p_values.csv')

        df_old_chi_squared.to_csv(f'{real_data_path}/{file_name}/results/haplotypes/old_chi_squared_statistic.csv')
        df_new_chi_squared.to_csv(f'{real_data_path}/{file_name}/results/haplotypes/new_chi_squared_statistic.csv')
        # df_ambiguities.to_csv(f'{real_data_path}/results/ambiguities')
        df_dof.to_csv(f'{real_data_path}/{file_name}/results/haplotypes/dof.csv')


def perform_tests_save_data(file_name, is_using_haplotypes, level_haplotype, race):
    # df_id_allele1_allele2_prob = pd.read_csv(
    #     f'{real_data_path}/levels/{level}/races/{race}/id_allele1_allele2_probability',
    #     dtype={'id': str, 'allele_1': int, 'allele_2': int, 'probability': float})
    id_to_index = {}
    allele_to_index = {}

    if not is_using_haplotypes:
        probabilities_path = f'{real_data_path}/{file_name}/levels/{level_haplotype}/races/{race}/id_allele1_allele2_probability'
    else:
        probabilities_path = f'{real_data_path}/{file_name}/haplotypes/races/{race}/id_allele1_allele2_probability'

    # first read all the rows and get indices of ids and alleles and amounts
    with open(probabilities_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # header
            if index == 0:
                continue
            if index % 10000 == 0:
                if not is_using_haplotypes:
                    print(f'row: {index}, file_name: {file_name}, level: {level_haplotype}, race: {race}, preprocessing')
                else:
                    print(f'row: {index}, level: {level_haplotype}, race: {race}, preprocessing')
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

    if is_using_haplotypes:
        print(f'race: {race}, not using haplotypes, alleles amount: {alleles_count}')

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
                if not is_using_haplotypes:
                    print(f'row: {index}, file_name: {file_name}, level: {level_haplotype}, race: {race}, After preprocessing')
                else:
                    print(f'row: {index}, level: {level_haplotype}, race: {race}, After preprocessing')
            lst = line.strip('\n').split(',')

            # id = lst[0]
            allele_1 = lst[1]
            allele_2 = lst[2]
            probability = float(lst[3])

            # id_index = id_to_index[id]

            allele_1_index = allele_to_index[allele_1]
            allele_2_index = allele_to_index[allele_2]

            # print(f'allele_1: {allele_1}, allele_2: {allele_2}')
            # print(f'allele_1_index: {allele_1_index}, allele_2_index: {allele_2_index}')
            allele_1_index, allele_2_index = min(allele_1_index, allele_2_index), max(allele_1_index, allele_2_index)
            # print(f'allele_1_index: {allele_1_index}, allele_2_index: {allele_2_index}')
            # print('----------------------------------------------')
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

    return old_p_value, new_p_value, old_chi_squared, new_chi_squared, dof


if __name__ == '__main__':
    file_names = constants.FILE_NAMES

    for file in file_names:
        # run for alleles
        get_data_for_chi_squared(file_name=file,
                                 is_using_haplotypes=0)
        # run for haplotypes
        # get_data_for_chi_squared(file_name=file,
        #                          is_using_haplotypes=1)
