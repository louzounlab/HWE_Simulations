import pandas as pd
import os

import utils


def create_alleles_indices_files():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    # level -> dict from alleles to indices
    dict_dicts_alleles_indices = {}

    columns = {'allele': pd.Series(dtype='str'),
               'index': pd.Series(dtype='float')}

    # level -> df [alleles, indices]
    dict_dfs_alleles_indices = {}

    with open("data/us_grma.haps.freqs", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            print(index)
            lst = line.split(',')
            # ['47919664', 'A*02:01~B*51:01~C*02:02~DQB1*03:01~DRB1*11:01', 'A*24:02~B*27:02~C*02:02~DQB1*03:01~DRB1*13:05', '2.9624149804921254e-08', '16\n']
            # '47919664'
            id = lst[0]
            # 'A*02:01~B*51:01~C*02:02~DQB1*03:01~DRB1*11:01'
            first_alleles = lst[1]
            # 'A*02:01~B*51:01~C*02:02~DQB1*03:01~DRB1*11:01'
            second_alleles = lst[2]
            # '2.16e-13'
            probability = lst[3]

            for level in levels_set:
                if level not in dict_dicts_alleles_indices:
                    # add new dictionary for the level
                    dict_dicts_alleles_indices[level] = {}
                    dict_dfs_alleles_indices[level] = pd.DataFrame(columns)
                dict_alleles_indices = dict_dicts_alleles_indices[level]
                df = dict_dfs_alleles_indices[level]

                level_position = utils.LEVELS_DICT[level]
                # [A*02:01, B*02:01, C*02:01]
                left_alleles = first_alleles.split('~')
                # [A*02:01, B*02:01, C*02:01]
                right_alleles = second_alleles.split('~')

                # [A*02, B*02, C*02]
                left_alleles = [sub.split(':')[0] for sub in left_alleles]
                right_alleles = [sub.split(':')[0] for sub in right_alleles]

                # 'A', '02'
                left_level, left_allele = left_alleles[level_position].split('*')
                # 'A', '03'
                right_level, right_allele = right_alleles[level_position].split('*')

                # put alleles
                if left_allele not in dict_alleles_indices:
                    dict_alleles_indices[left_allele] = len(dict_alleles_indices)
                    df.loc[len(dict_alleles_indices) - 1] = [left_allele, len(dict_alleles_indices) - 1]
                if right_allele not in dict_alleles_indices:
                    dict_alleles_indices[right_allele] = len(dict_alleles_indices)
                    df.loc[len(dict_alleles_indices) - 1] = [right_allele, len(dict_alleles_indices) - 1]
    for level in levels_set:
        dict_dfs_alleles_indices[level].to_csv(f'data/real_data/levels/{level}/alleles_indices', index=False)


# for every level and race creates an id, amount file
def create_amounts_files():
    races_set = utils.RACES_SET

    race_to_id_sum = {}

    current_race = ''
    current_id = ''
    current_sum = 0.0

    races_df = pd.read_csv('data/id_race_don_us.csv',
                           dtype={'id': str, 'rcae': str})
    # id -> race
    races_dict = {}
    for i in range(len(races_df)):
        if i % 500 == 0:
            print(f'Creating races dictionary: {i}')
        id = races_df.loc[i, 'id']
        race = races_df.loc[i, 'rcae']
        if id not in races_dict:
            races_dict[id] = race

    with open("data/us_grma.haps.freqs", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            print(index)
            lst = line.split(',')

            id_row = lst[0]
            observed_row = float(lst[3])
            person_index_row = int(lst[4])

            if person_index_row == 0 and current_id:
                # save the last person
                if current_race not in race_to_id_sum:
                    race_to_id_sum[current_race] = []
                # race_to_id_amount[current_race][current_id] = current_amount
                race_to_id_sum[current_race].append({'id': current_id, 'sum': current_sum})
                current_sum = 0.0

            current_race = races_dict[id_row]
            current_id = id_row
            current_sum += observed_row

    # save also the last person
    if current_race not in race_to_id_sum:
        race_to_id_sum[current_race] = []
    race_to_id_sum[current_race].append({'id': current_id, 'sum': current_sum})

    # save the files
    for race in races_set:
        list_of_dicts = race_to_id_sum[race]
        df = pd.DataFrame.from_dict(list_of_dicts)
        df.to_csv(f'data/real_data/races/{race}/id_sum', index=False)


# every person appears in subsequent rows with an index indicating how many times he appears.
# for every level and race we create a probabilities file that contains ID,ALLELE_1,ALLELE_2,O_K(I,J)
def create_id_allele1_allele2_probability_files_for_each_level(level):
    races_set = utils.RACES_SET

    # dictionary (race -> dict of probabilities [ID,ALLELE_1,ALLELE_2,O_K(I,J)] )
    race_to_list_of_dicts = {}

    # dictionary (allele -> index of allele)
    allele_index_dict = {}
    dtypes = {'allele': str, 'index': int}
    df = pd.read_csv(f'data/real_data/levels/{level}/alleles_indices',
                     dtype=dtypes)
    # allele: [0], allele: [1],...
    for i in range(len(df)):
        print(f'Creating alleles indices dictionary: {i}')
        allele = df.loc[i, 'allele']
        index = df.loc[i, 'index']
        if allele not in allele_index_dict:
            allele_index_dict[allele] = index
    # allele: 0, allele: 1,...

    races_df = pd.read_csv('data/id_race_don_us.csv',
                           dtype={'id': str, 'rcae': str})
    # id -> race
    races_dict = {}
    for i in range(len(races_df)):
        if i % 500 == 0:
            print(f'Creating races dictionary: {i}')
        id = races_df.loc[i, 'id']
        race = races_df.loc[i, 'rcae']
        if id not in races_dict:
            races_dict[id] = race


    # race -> dictionary (id -> amount)
    race_to_id_sum = {}
    for race in races_set:
        dtypes = {'id': str, 'sum': float}
        df = pd.read_csv(f'data/real_data/races/{race}/id_sum',
                         dtype=dtypes)
        my_dict = {}
        for i in range(len(df)):
            print(f'Creating id amount dictionary: {i}')
            id = df.loc[i, 'id']
            sum = df.loc[i, 'sum']
            my_dict[id] = sum
        race_to_id_sum[race] = my_dict

    race_to_id_allele1_allele2_to_probability = {}

    with open("data/us_grma.haps.freqs", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 1000 == 0:
                print(f'index: {index}, level: {level}')
            lst = line.split(',')
            # ['47919664', 'A*02:01~B*51:01~C*02:02~DQB1*03:01~DRB1*11:01', 'A*24:02~B*27:02~C*02:02~DQB1*03:01~DRB1*13:05', '2.9624149804921254e-08', '16\n']
            # '47919664'
            id = lst[0]
            # 'A*02:01~B*51:01~C*02:02~DQB1*03:01~DRB1*11:01'
            first_alleles = lst[1]
            # 'A*02:01~B*51:01~C*02:02~DQB1*03:01~DRB1*11:01'
            second_alleles = lst[2]
            # '2.16e-13'
            probability = lst[3]

            person_index = lst[4]

            # [A*02:01, B*02:01, C*02:01]
            left_alleles = first_alleles.split('~')
            # [A*02:01, B*02:01, C*02:01]
            right_alleles = second_alleles.split('~')

            # [A*02, B*02, C*02]
            left_alleles = [sub.split(':')[0] for sub in left_alleles]
            right_alleles = [sub.split(':')[0] for sub in right_alleles]

            race = races_dict[id]
            if race == 'UNK':
                continue
            probability = float(probability) / race_to_id_sum[race][id]

            level_position = utils.LEVELS_DICT[level]
            # 'A', '02'
            left_level, left_allele = left_alleles[level_position].split('*')
            # 'A', '03'
            right_level, right_allele = right_alleles[level_position].split('*')

            # # dataframe: allele -> index of allele
            left_allele_index = allele_index_dict[left_allele]
            right_allele_index = allele_index_dict[right_allele]

            # # level -> dictionary (race -> df of probabilities [ID,ALLELE_1,ALLELE_2,O_K(I,J)] )
            # dict_levels_races_probabilities = {}
            if race not in race_to_list_of_dicts:
                race_to_list_of_dicts[race] = []
            list_of_dicts = race_to_list_of_dicts[race]
            # observed_df.loc[len(observed_df)] = [id, left_allele_index, right_allele_index, float(observed)]

            if race not in race_to_id_allele1_allele2_to_probability:
                race_to_id_allele1_allele2_to_probability[race] = {}

            left, right = min(left_allele_index, right_allele_index), max(left_allele_index, right_allele_index)
            id_allele1_allele2 = id + str(left) + str(right)

            # if its the first time we observe this person with those alleles
            if id_allele1_allele2 not in race_to_id_allele1_allele2_to_probability[race]:
                race_to_id_allele1_allele2_to_probability[race][id_allele1_allele2] = probability
                list_of_dicts.append({'id': id, 'allele_1': left_allele_index, 'allele_2': right_allele_index,
                                      'probability': float(probability)})
            else:
                # we already observed this person with those alleles
                list_of_dicts[len(list_of_dicts) - 1]['probability'] += probability

    for race in races_set:
        list_of_dicts = race_to_list_of_dicts[race]
        df = pd.DataFrame.from_dict(list_of_dicts)
        df.to_csv(f'data/real_data/levels/{level}/races/{race}/id_allele1_allele2_probability', index=False)


# in B level there are races with one less
def check_alleles_indices(level):
    races_set = utils.RACES_SET
    for race in races_set:
        alleles_set = set()
        dtypes = {'id': str, 'allele_1': int, 'allele_2': int, 'probability': float}
        df = pd.read_csv(f'data/real_data/levels/{level}/races/{race}/id_allele1_allele2_probability',
                         dtype=dtypes)
        for index, row in df.iterrows():
            left_allele = row['allele_1']
            right_allele = row['allele_2']
            if left_allele not in alleles_set:
                alleles_set.add(left_allele)
            if right_allele not in alleles_set:
                alleles_set.add(right_allele)
        print(f'amount of alleles in level {level}, race {race} is: {len(alleles_set)}')


# for each level and race we create a csv file 'uncertainty' with id,uncertainty
def create_uncertainty_file_for_each_level_race():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    for level in levels_set:
        for race in races_set:
            print(f'level: {level}, race: {race}')
            # read id_allele1_allele2_probability file
            # for every different ID calculate 1-sum_{i,j}p_k(i,j)^2
            df = pd.read_csv(f'data/real_data/levels/{level}/races/{race}/id_allele1_allele2_probability',
                             dtype={'id': str, 'allele_1': int, 'allele_2': int, 'probability': float})
            uncertainty_list_of_dicts = []
            last_id = ''
            last_uncertainty = 1.0
            for index, row in df.iterrows():
                current_id = row['id']
                current_probability = row['probability']

                # finished going over all the rows of the last id
                if (current_id != last_id) and last_id:
                    uncertainty_list_of_dicts.append({'id': last_id, 'uncertainty': last_uncertainty})
                    last_uncertainty = 1.0

                last_id = current_id
                last_uncertainty -= (current_probability ** 2)

            uncertainty_df = pd.DataFrame.from_dict(uncertainty_list_of_dicts)
            uncertainty_df.to_csv(f'data/real_data/levels/{level}/races/{race}/uncertainty', index=False)


# create i_j_probability files
def create_probabilities_file_for_each_level_race():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    amounts_dict_dict = utils.AMOUNTS_DICT_DICTS

    for level in levels_set:
        for race in races_set:
            print(f'level: {level}, race: {race}')
            # i_j -> {'i':, 'j':, 'probability': }
            i_j_to_dict = {}
            allele_to_index = {}
            df = pd.read_csv(f'data/real_data/levels/{level}/races/{race}/id_allele1_allele2_probability',
                             dtype={'id': str, 'allele_1': int, 'allele_2': int, 'probability': float})
            for index, row in df.iterrows():
                left_allele = row['allele_1']
                print(left_allele)
                right_allele = row['allele_2']
                probability = row['probability']

                if left_allele not in allele_to_index:
                    allele_to_index[left_allele] = len(allele_to_index)
                if right_allele not in allele_to_index:
                    allele_to_index[right_allele] = len(allele_to_index)

                i = allele_to_index[left_allele]
                j = allele_to_index[right_allele]
                i, j = min(i, j), max(i, j)

                i_j = str(i) + '_' + str(j)
                if i_j not in i_j_to_dict:
                    i_j_to_dict[i_j] = {'i': i, 'j': j, 'probability': probability}
                else:
                    i_j_to_dict[i_j]['probability'] += probability
            # we have a dictionary: i_j -> {'i':, 'j':, 'probability': }
            # translate it to a list of dictionaries:
            # but first normalize p(i,j) = 1/N * sum_k p_k(i,j)
            for key in i_j_to_dict:
                i_j_to_dict[key]['probability'] /= amounts_dict_dict[level][race]
            list_of_dicts = list(i_j_to_dict.values())

            df_probabilities = pd.DataFrame.from_dict(list_of_dicts)
            df_probabilities.to_csv(f'data/real_data/levels/{level}/races/{race}/i_j_probability', index=False)


# for every level and race create a csv file of i, probability
def create_alleles_probability_files():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    for level in levels_set:
        for race in races_set:
            print(f'level: {level}, race: {race}')
            i_to_i_probability = {}

            df = pd.read_csv(f'data/real_data/levels/{level}/races/{race}/i_j_probability',
                             dtype={'i': int, 'j': int, 'probability': float})
            for index, row in df.iterrows():
                i = int(row['i'])
                j = int(row['j'])
                probability = row['probability']

                if i not in i_to_i_probability:
                    i_to_i_probability[i] = {'i': i, 'probability': probability / 2}
                else:
                    i_to_i_probability[i]['probability'] += probability / 2

                if j not in i_to_i_probability:
                    i_to_i_probability[j] = {'i': j, 'probability': probability / 2}
                # if j is already in the dictionary and is not i (so we don't count probability twice)
                else:
                    i_to_i_probability[j]['probability'] += probability / 2

            list_of_dicts = list(i_to_i_probability.values())

            df_probabilities = pd.DataFrame.from_dict(list_of_dicts)
            df_probabilities.to_csv(f'data/real_data/levels/{level}/races/{race}/i_probability', index=False)


