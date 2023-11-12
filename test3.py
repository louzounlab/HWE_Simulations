import pandas as pd
import numpy as np
import csv

# level = 'A'
# race = 'AFA'
#
# path = f'data/real_data_test/levels/{level}/races/{race}'
#
# ambiguity_df = pd.read_csv(f'{path}/uncertainty',
#                            dtype={'id': str, 'uncertainty': float})
# i_j_probability_df = pd.read_csv(f'{path}/id_allele1_allele2_probability',
#                                  dtype={'id': str, 'allele_1': int, 'allele_2': int, 'probability': float})
# id_sum_df = pd.read_csv('data/real_data_test/races/AFA/id_sum',
#                         dtype={'id': str, 'sum': float})
# print(ambiguity_df.shape[0])
# print(i_j_probability_df.shape[0])
# print(id_sum_df.shape[0])
#
# i_j_probability_B_df = pd.read_csv(f'data/real_data_test/levels/B/races/{race}/id_allele1_allele2_probability',
#                                  dtype={'id': str, 'allele_1': int, 'allele_2': int, 'probability': float})
# print(i_j_probability_B_df.shape[0])

# first read all the rows and get indices of ids and alleles and amounts
# with open('profiles1.csv', encoding="utf8") as infile:
#     reader = csv.reader(infile)
#     lst = next(reader)
#
#     name_col = lst[0]
#     age_col = lst[1]
#     country_col = lst[2]
#
# df = pd.read_csv('profiles1.csv',
#                  dtype={name_col: str, age_col: int, country_col: str})
# id_col = df.columns[0]
# df = df.sort_values(by=id_col)
# print(df.loc[0])
# print(df.columns[0])

level = 'A'
race = 'HIS'
from models import chi_squared
real_data_path = 'data/real_data_test'
probabilities_path = f'{real_data_path}/levels/{level}/races/{race}/id_allele1_allele2_probability'

id_to_index = {}
allele_to_index = {}
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
        allele_1 = int(lst[1])
        allele_2 = int(lst[2])
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

        id = lst[0]
        allele_1 = int(lst[1])
        allele_2 = int(lst[2])
        probability = float(lst[3])

        id_index = id_to_index[id]

        allele_1_index = allele_to_index[allele_1]
        allele_2_index = allele_to_index[allele_2]

        alleles_probabilities[allele_1_index] += 0.5 * probability
        alleles_probabilities[allele_2_index] += 0.5 * probability

        observed_probabilities[allele_1_index, allele_2_index] += probability

        correction[allele_1_index, allele_2_index] += (probability ** 2)

alleles_probabilities /= population_amount

for i in range(alleles_count):
    for j in range(alleles_count):
        observed_probabilities[i, j] /= population_amount
        if observed_probabilities[i, j] == 0:
            correction[i, j] = 1.0
        else:
            correction[i, j] /= (population_amount * observed_probabilities[i, j])

min_i = 0
min_j = 0
min_corr = correction[0, 0]

for i in range(alleles_count):
    for j in range(i, alleles_count):
        expected = population_amount * alleles_probabilities[i] * alleles_probabilities[j]
        if i != j:
            expected *= 2
        if expected >= 2.0 and correction[i, j] < min_corr:
            min_corr = correction[i, j]
            min_i = i
            min_j = j

print(min_i)
print(min_j)

expected = population_amount * alleles_probabilities[min_i] * alleles_probabilities[min_j]
observed = population_amount * observed_probabilities[min_i, min_j]

print(f'corr: {min_corr}, expected: {expected}, observed: {observed}')