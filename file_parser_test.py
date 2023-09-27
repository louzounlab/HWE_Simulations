import pandas as pd
import utils

real_data_path = 'data/real_data_test'


# helper function. from '05:09N' to '0509'
def allele_format_to_number_string(s):
    s_new = ''
    for char in s:
        if char.isdigit():
            s_new += char
        # if char.isalpha():
        #     s_new += str(ord(char))
    return s_new


# for every level and race creates an id, sum file
def create_sums_files():
    races_set = utils.RACES_SET

    race_to_id_sum = {}

    current_race = ''
    current_id = ''
    current_sum = 0.0

    races_df = pd.read_csv('data/id_race_don_us.csv',
                           dtype={'id': str, 'rcae': str})

    # id -> race
    races_dict = {}
    for index, row in races_df.iterrows():
        id_row = row['id']
        race_row = row['rcae']
        if id_row not in races_dict:
            races_dict[id_row] = race_row

    with open("data/us_grma.haps.freqs", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 500 == 0:
                print(index)
            lst = line.split(',')

            id_row = lst[0]
            observed_row = float(lst[3])
            person_index_row = int(lst[4])

            # if we got to a new id
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
        df.to_csv(f'{real_data_path}/races/{race}/id_sum', index=False)


# every person appears in subsequent rows with an index indicating how many times he appears.
# for every level and race we create a probabilities file that contains ID,ALLELE_1,ALLELE_2,O_K(I,J)
def create_id_allele1_allele2_probability_files_for_each_level(level):
    # dictionary (race -> dict of probabilities [ID,ALLELE_1,ALLELE_2,O_K(I,J)] )
    race_to_list_of_dicts = {}

    races_df = pd.read_csv('data/id_race_don_us.csv',
                           dtype={'id': str, 'rcae': str})
    races_set = utils.RACES_SET
    # id -> race
    races_dict = {}
    for index, row in races_df.iterrows():
        if index % 500 == 0:
            print(f'Creating races dictionary: {index}')
        id = row['id']
        race = row['rcae']
        if id not in races_dict:
            races_dict[id] = race

    # race -> dictionary (id -> sum)
    race_to_id_to_sum = {}
    for race in races_set:
        dtypes = {'id': str, 'sum': float}
        df = pd.read_csv(f'{real_data_path}/races/{race}/id_sum',
                         dtype=dtypes)
        my_dict = {}
        for index, row in df.iterrows():
            print(f'Creating id amount dictionary: {index}')
            id_row = row['id']
            sum_row = row['sum']
            my_dict[id_row] = sum_row
        race_to_id_to_sum[race] = my_dict

    race_to_id_allele1_allele2_to_probability = {}

    # race -> id_allele1_allele2 -> dict {id: allele_1: allele_2: probability:}
    race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict = {}

    with open("data/us_grma.haps.freqs", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 1000 == 0:
                print(f'index: {index}, level: {level}')
            lst = line.split(',')
            # ['47919664', 'A*02:01~B*51:01~C*02:02~DQB1*03:01~DRB1*11:01', 'A*24:02~B*27:02~C*02:02~DQB1*03:01~DRB1*13:05', '2.9624149804921254e-08', '16\n']
            # '47919664'
            id_row = lst[0]
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

            level_position = utils.LEVELS_DICT[level]

            # 'A', '02:01N'
            left_level, left_allele = left_alleles[level_position].split('*')
            right_level, right_allele = right_alleles[level_position].split('*')

            # '0201'
            left_allele = allele_format_to_number_string(left_allele)
            right_allele = allele_format_to_number_string(right_allele)

            race = races_dict[id_row]
            if race == 'UNK':
                continue
            probability = float(probability) / race_to_id_to_sum[race][id_row]

            if race not in race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict:
                race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race] = {}

            left, right = min(int(left_allele), int(right_allele)), max(int(left_allele), int(right_allele))
            id_allele1_allele2 = id_row + '_' + str(left) + '_' + str(right)

            # if its the first time we observe this person with those alleles
            if id_allele1_allele2 not in race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race]:
                race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race][id_allele1_allele2] = {
                    'id': id_row,
                    'allele_1': left,
                    'allele_2': right,
                    'probability': float(probability)
                }
            else:
                # we already observed this person with those alleles
                race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race][id_allele1_allele2][
                    'probability'] += probability

    for race in races_set:
        # id_allele1_allele2 -> {id, allele_1, allele_2, probability}
        my_dict = race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race]
        # create a dataframe from a list of dicts
        df = pd.DataFrame.from_dict(list(my_dict.values()))
        df.to_csv(f'{real_data_path}/levels/{level}/races/{race}/id_allele1_allele2_probability', index=False)


#### FINE UNTIL HERE


# for each level and race we create a csv file 'uncertainty' with id,uncertainty
def create_uncertainty_file_for_each_level_race():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    for level in levels_set:
        for race in races_set:
            print(f'level: {level}, race: {race}')
            # read id_allele1_allele2_probability file
            # for every different ID calculate 1-sum_{i,j} p_k(i,j)^2
            df = pd.read_csv(f'{real_data_path}/levels/{level}/races/{race}/id_allele1_allele2_probability',
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

            # adding uncertainty for the last person
            uncertainty_list_of_dicts.append({'id': last_id, 'uncertainty': last_uncertainty})

            # saving csv
            uncertainty_df = pd.DataFrame.from_dict(uncertainty_list_of_dicts)
            uncertainty_df.to_csv(f'{real_data_path}/levels/{level}/races/{race}/uncertainty', index=False)


def check_amounts():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    for level in levels_set:
        for race in races_set:
            df = pd.read_csv(f'{real_data_path}/levels/{level}/races/{race}/uncertainty')
            print(f'level: {level}, race: {race}, population: {df.shape[0]}')


# create i_j_probability files
def create_probabilities_file_for_each_level_race():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    # level -> race -> sums
    amounts_dict_dict = utils.AMOUNTS_DICT_DICTS

    for level in levels_set:
        for race in races_set:
            print(f'level: {level}, race: {race}')
            # i_j -> {'i':, 'j':, 'probability': }
            i_j_to_dict = {}
            allele_to_index = {}
            df = pd.read_csv(f'{real_data_path}/levels/{level}/races/{race}/id_allele1_allele2_probability',
                             dtype={'id': str, 'allele_1': int, 'allele_2': int, 'probability': float})
            for index, row in df.iterrows():
                left_allele = row['allele_1']
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

            for key in i_j_to_dict:
                i_j_to_dict[key]['probability'] /= amounts_dict_dict[level][race]
            list_of_dicts = list(i_j_to_dict.values())

            df_probabilities = pd.DataFrame.from_dict(list_of_dicts)
            df_probabilities.to_csv(f'{real_data_path}/levels/{level}/races/{race}/i_j_probability', index=False)


# for every level and race create a csv file of i, probability
def create_alleles_probability_files():
    levels_set = utils.LEVELS_SET
    races_set = utils.RACES_SET

    for level in levels_set:
        for race in races_set:
            print(f'level: {level}, race: {race}')
            i_to_i_probability = {}

            df = pd.read_csv(f'{real_data_path}/levels/{level}/races/{race}/i_j_probability',
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
            df_probabilities.to_csv(f'{real_data_path}/levels/{level}/races/{race}/i_probability', index=False)


create_alleles_probability_files()
