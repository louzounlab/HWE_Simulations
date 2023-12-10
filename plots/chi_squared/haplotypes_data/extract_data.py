import pandas as pd
import constants

real_data_path = '../../../data/haplotypes_data'

# THIS FILE EXTRACT THE HAPLOTYPES DATA AND SAVE IN THE DIRECTORIES


# for every race creates an id, sum file
def create_sums_files(file_name):
    races_list = constants.FILE_NAME_TO_POPULATIONS[file_name]
    race_to_id_sum = {}

    current_race = ''
    current_id = ''
    current_sum = 0.0

    # access the races
    # id -> race
    races_dict = {}
    with open(f'{real_data_path}/{file_name}_for_impute_nopyard.csv', encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 500 == 0:
                print(f'Reading races and ids. index: {index}')
            # ['D009999', 'A*11:01~C*14:02~B*51:01~DRB1*15:01~DQB1*06:02+A*24:07~C*04:01~B*35:05~DRB1*09:01~DQB1*03:03', '1', 'FILII']
            lst = line.strip('\n').split(',')
            id_row = lst[0]
            race_row = lst[3]

            if id_row not in races_dict:
                races_dict[id_row] = race_row

    # calculating the sums for each id
    with open(f'{real_data_path}/{file_name}_nopyard.freqs', encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 500 == 0:
                print(index)
            # ['D009998', 'A*11:01~B*48:01~C*08:01~DQB1*03:01~DRB1*11:01+A*24:07~B*18:01~C*07:02~DQB1*05:01~DRB1*15:02', '2.100112733255187e-08', '1']
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

    # save the files (for each single population)
    list_of_dicts_all_races = []
    for race in races_list:
        if race == 'ALL':
            continue
        list_of_dicts = race_to_id_sum[race]
        list_of_dicts_all_races += list_of_dicts
        df = pd.DataFrame.from_dict(list_of_dicts)
        df.to_csv(f'{real_data_path}/{file_name}/races/{race}/id_sum', index=False)

    # save the files for all population
    df = pd.DataFrame.from_dict(list_of_dicts_all_races)
    df.to_csv(f'{real_data_path}/{file_name}/races/ALL/id_sum', index=False)


# every person appears in subsequent rows with an index indicating how many times he appears.
# for every level and race we create a probabilities file that contains ID,ALLELE_1,ALLELE_2,O_K(I,J)
def create_id_allele1_allele2_probability_files(file_name, is_using_haplotypes, level_haplotype):
    # file_name: Name of the file (without .csv)
    # is_using_haplotypes: whether we use the alleles or the haplotypes
    # level_haplotype: either the level we use, or None if we use haplotypes

    # dictionary (race -> dict of probabilities [ID,ALLELE_1,ALLELE_2,O_K(I,J)] )
    # race_to_list_of_dicts = {}

    # races_df = pd.read_csv('data/id_race_don_us.csv',
    #                        dtype={'id': str, 'rcae': str})
    races_list = constants.FILE_NAME_TO_POPULATIONS[file_name]
    # access the races
    # id -> race
    races_dict = {}
    with open(f'{real_data_path}/{file_name}_for_impute_nopyard.csv', encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 500 == 0:
                print(f'Reading races and ids. index: {index}')
            # ['D009999', 'A*11:01~C*14:02~B*51:01~DRB1*15:01~DQB1*06:02+A*24:07~C*04:01~B*35:05~DRB1*09:01~DQB1*03:03', '1', 'FILII']
            lst = line.strip('\n').split(',')
            id_row = lst[0]
            race_row = lst[3]

            if id_row not in races_dict:
                races_dict[id_row] = race_row

    # race -> dictionary (id -> sum)
    print(f'Creating id amount dictionary:')
    race_to_id_to_sum = {}
    for race in races_list:
        # dtypes = {'id': str, 'sum': float}
        # df = pd.read_csv(f'{real_data_path}/races/{race}/id_sum',
        #                  dtype=dtypes)
        my_dict = {}
        with open(f'{real_data_path}/{file_name}/races/{race}/id_sum', encoding="utf8") as infile:
            for index, line in enumerate(infile):
                if index == 0:
                    continue
                lst = line.strip('\n').split(',')
                id_row = lst[0]
                sum_row = lst[1]
                my_dict[id_row] = float(sum_row)
        race_to_id_to_sum[race] = my_dict

    # FINE UNTIL HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # race -> id_allele1_allele2 -> dict {id: allele_1: allele_2: probability:}
    race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict = {}

    with open(f'{real_data_path}/{file_name}_nopyard.freqs', encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 1000 == 0:
                if not is_using_haplotypes:
                    print(f'index: {index}, file_name: {file_name}, level: {level_haplotype}')
                else:
                    print(f'index: {index}, file_name: {file_name}, haplotypes')
            lst = line.split(',')
            # '47919664'
            id_row = lst[0]

            # '2.16e-13'
            probability = float(lst[2])

            # getting the alleles or haplotypes
            if not is_using_haplotypes:
                level_position = constants.LEVELS_DICT[level_haplotype]
                # ['A*24:07~B*48:01~C*08:01~DQB1*03:01~DRB1*11:01', 'A*11:01~B*18:01~C*07:02~DQB1*05:01~DRB1*15:02']
                alleles_splitted = lst[1].split('+')
                left_allele = alleles_splitted[0].split('~')[level_position]
                right_allele = alleles_splitted[1].split('~')[level_position]
            else:
                alleles_splitted = lst[1].split('+')
                left_allele = alleles_splitted[0]
                right_allele = alleles_splitted[1]

            race = races_dict[id_row]

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

    my_dict_for_all_population = {}
    # saving for single populations
    for race in races_list:
        if race == 'ALL':
            continue
        # id_allele1_allele2 -> {id, allele_1, allele_2, probability}
        my_dict = race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race]
        # merging current dictionary to the all population dictionary
        my_dict_for_all_population.update(my_dict)
        # create a dataframe from a list of dicts
        df = pd.DataFrame.from_dict(list(my_dict.values()))
        if not is_using_haplotypes:
            df.to_csv(f'{real_data_path}/{file_name}/levels/{level_haplotype}/races/{race}/id_allele1_allele2_probability', index=False)
        else:
            df.to_csv(f'{real_data_path}/{file_name}/haplotypes/races/{race}/id_allele1_allele2_probability', index=False)

    # saving for all the population
    df = pd.DataFrame.from_dict(list(my_dict_for_all_population.values()))
    if not is_using_haplotypes:
        df.to_csv(f'{real_data_path}/{file_name}/levels/{level_haplotype}/races/ALL/id_allele1_allele2_probability', index=False)
    else:
        df.to_csv(f'{real_data_path}/{file_name}/haplotypes/races/ALL/id_allele1_allele2_probability', index=False)


if __name__ == '__main__':
    level_list = constants.LEVELS_LIST
    file_names = constants.FILE_NAMES

    for file in file_names:
        create_sums_files(file_name=file)

    for file in file_names:
        for level in level_list:
            create_id_allele1_allele2_probability_files(file_name=file,
                                                        is_using_haplotypes=0,
                                                        level_haplotype=level)
        # run also for haplotypes
        create_id_allele1_allele2_probability_files(file_name=file,
                                                    is_using_haplotypes=1,
                                                    level_haplotype=None)
