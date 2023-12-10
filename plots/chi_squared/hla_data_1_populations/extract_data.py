import pandas as pd
import utils

real_data_path = '../../../data/hla_data_1_populations'
us_grma = '../../../data/us_grma_pyard.freqs'
# real_data_path = 'data/real_data_test'


# # helper function. from 'C*05:09N' to 'C*05:09'
# def allele_format_to_number_string(s):
#     # take the prefix 'C*'
#     s_new = s[:2]
#     for char in s[2:]:
#         if not char.isalpha():
#             s_new += char
#         # if char.isdigit():
#         #     s_new += char
#     return s_new


# for every level and race creates an id, sum file
def create_sums_files():
    races_set = utils.RACES_1_LIST

    race_to_id_sum = {}

    # current_race = ''
    current_id = ''
    current_sum = 0.0

    # races_df = pd.read_csv('data/id_race_don_us.csv',
    #                        dtype={'id': str, 'rcae': str})

    # id -> race
    races_dict = {}
    with open('../../../data/id_race_don_us.csv', encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 500 == 0:
                print(f'Reading races and ids. index: {index}')
            if index == 0:
                continue
            lst = line.strip('\n').split(',')

            id_row = lst[0]
            race_row = lst[1]
            if race_row != 'UNK':
                race_row = 'ALL'
            if id_row not in races_dict:
                races_dict[id_row] = race_row

    with open(us_grma, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 500 == 0:
                print(index)
            lst = line.strip('\n').split(',')

            id_row = lst[0]
            observed_row = float(lst[2])
            person_index_row = int(lst[3])

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


# create_sums_files()


# every person appears in subsequent rows with an index indicating how many times he appears.
# for every level and race we create a probabilities file that contains ID,ALLELE_1,ALLELE_2,O_K(I,J)
def create_id_allele1_allele2_probability_files_for_each_level(level_):
    # dictionary (race -> dict of probabilities [ID,ALLELE_1,ALLELE_2,O_K(I,J)] )
    # race_to_list_of_dicts = {}

    # races_df = pd.read_csv('data/id_race_don_us.csv',
    #                        dtype={'id': str, 'rcae': str})
    races_set = utils.RACES_1_LIST
    # id -> race
    races_dict = {}
    print(f'Creating races dictionary:')
    with open("../../../data/id_race_don_us.csv", encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index == 0:
                continue
            lst = line.strip('\n').split(',')

            id_row = lst[0]
            race_row = lst[1]
            if id_row not in races_dict:
                races_dict[id_row] = race_row

    # race -> dictionary (id -> sum)
    print(f'Creating id amount dictionary:')
    race_to_id_to_sum = {}
    for race in races_set:
        # dtypes = {'id': str, 'sum': float}
        # df = pd.read_csv(f'{real_data_path}/races/{race}/id_sum',
        #                  dtype=dtypes)
        my_dict = {}
        with open(f'{real_data_path}/races/{race}/id_sum', encoding="utf8") as infile:
            for index, line in enumerate(infile):
                if index == 0:
                    continue
                lst = line.strip('\n').split(',')
                id_row = lst[0]
                sum_row = lst[1]
                my_dict[id_row] = float(sum_row)
        race_to_id_to_sum[race] = my_dict

    # FIXED UNTIL HERE
    # race -> id_allele1_allele2 -> dict {id: allele_1: allele_2: probability:}
    race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict = {}

    with open(us_grma, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 1000 == 0:
                print(f'index: {index}, level_: {level_}')
            lst = line.split(',')
            # '47919664'
            id_row = lst[0]

            # '2.16e-13'
            probability = float(lst[2])

            level_position = utils.LEVELS_DICT[level_]

            # ['A*11:01+A*32:01, B*18:01+B*44:03, C*04:01+C*07:01, DQB1*02:01+DQB1*03:01, DRB1*07:01+DRB1*11:43']
            alleles_full_string = lst[1].split('^')

            left_allele = alleles_full_string[level_position].split('+')[0]
            right_allele = alleles_full_string[level_position].split('+')[1]

            race = races_dict[id_row]
            if race != 'UNK':
                race = 'ALL'
            if race not in races_set:
                continue
            probability = float(probability) / race_to_id_to_sum[race][id_row]

            if race not in race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict:
                race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race] = {}

            left, right = min(left_allele, right_allele), max(left_allele, right_allele)
            # left, right = min(int(left_allele), int(right_allele)), max(int(left_allele), int(right_allele))
            id_allele1_allele2 = id_row + '_' + str(left) + '_' + str(right)

            # if it's the first time we observe this person with those alleles
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
        df.to_csv(f'{real_data_path}/levels/{level_}/races/{race}/id_allele1_allele2_probability', index=False)


#### FINE UNTIL HERE
levels_set = utils.LEVELS_LIST
for level in levels_set:
    create_id_allele1_allele2_probability_files_for_each_level(level)