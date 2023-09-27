import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import random

import utils
import utils_with_certainty


# def binary_search(cdf, target: float = 0):
#     start = 0
#     end = len(cdf) - 1
#     while end > start + 1:
#         mid = (end + start) // 2
#         # print(f'start: {start}. mid: {mid}. end: {end}')
#         if cdf[mid] < target:
#             start = mid
#         elif (cdf[mid] >= target) and (mid > 0) and (cdf[mid - 1] >= target):
#             end = mid
#         elif cdf[mid] >= target:
#             return mid
#     # now start and end pointing to the elements closest to 0
#     # pick the index of the closer one
#     return end
#
#
# test_space = [([0.3, 0.4, 0.5], 0.49, 2),
#               ([0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.7, 1.0, 1.0], 0.4001, 1),
#               ([0.3, 0.4, 0.5], 0.5, 2),
#               ([0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.7, 1.0, 1.0], 0.5, 1),
#               ([0.7, 0.7, 0.7, 0.8, 0.81, 0.82, 0.83, 0.83, 1.0], 0.83, 6),
#               ([0.7, 0.7, 0.7, 0.8, 0.81, 0.82, 0.83, 0.83, 1.0], 0.820000005, 6),
#               ([0.3, 0.3, 0.5, 0.6, 0.6, 0.6], 0.6, 3),
#               ([0.3, 0.3, 0.5, 0.6, 0.6, 0.6], 0.59, 3)
#               ]
# for i, test in enumerate(test_space):
#     index = binary_search(test[0], target=test[1])
#     if index != test[2]:
#         print('nope')

# data = pd.read_csv('data/means_alpha=00', index_col=0)
# plt.figure(figsize=(12, 10))
# sns.heatmap(data, annot=True, cmap=plt.cm.CMRmap_r)
# plt.xlabel('Population')
# plt.ylabel('Alleles')
# plt.title('Mean of % positive results')
# plt.show()

# data = pd.read_csv('data/stds_alpha=00', index_col=0)
# plt.figure(figsize=(12, 10))
# sns.heatmap(data, annot=True, cmap=plt.cm.CMRmap_r)
# plt.xlabel('Population')
# plt.ylabel('Alleles')
# plt.title('Stds of % positive results')
# plt.show()

# data = pd.read_csv('data/positives_for_alpha_values_alleles=100_population=10000000', index_col=0)
# alleles_count_ = 100
# population_amount_ = 10000000
# interval_for_alpha = 0.125
#
# plot_means = data.iloc[:, 0]
# plot_stds = data.iloc[:, 1]
# print(plot_means)
# print(type(plot_means))
#
# alpha_values = np.arange(start=0.0, stop=1.0 + interval_for_alpha, step=interval_for_alpha)
# plt.errorbar(alpha_values, plot_means, plot_stds)
# plt.xlabel('Alpha values')
# plt.ylabel('% Positives')
# plt.title(f'Alleles: {alleles_count_}. Population: {population_amount_}')
# plt.show()

def search(cdf, target: float = 0):
    for i in range(len(cdf)):
        if cdf[i] >= target:
            return i


def binary_search(cdf, target: float = 0):
    start = 0
    end = len(cdf)
    while start < end:
        mid = (end + start) // 2
        # print(f'start: {start}. mid: {mid}. end: {end}')
        if cdf[mid] < target:
            start = mid + 1
        else:
            end = mid
    # now start and end pointing to the elements closest to 0
    # pick the index of the closer one
    return start


# dtypes = {'allele': str, 'index': int}
# df = pd.read_csv('data/real_data/levels/A/alleles_indices',
#                  dtype=dtypes)
# print(df.set_index('allele').T.to_dict('list'))

# races_df = pd.read_csv('data/id_race_don_us.csv',
#                        dtype={'id': str, 'rcae': str})
# print(races_df.head())
# races_dict = {}
# for i in range(len(races_df)):
#     print(i)
#     id = races_df.loc[i, 'id']
#     race = races_df.loc[i, 'rcae']
#     if id not in races_dict:
#         races_dict[id] = race
# print(races_dict)

# df = pd.read_csv('data/id_race_don_us.csv')
# ids_amount = len(df)
# print(f'amount of ids: {ids_amount}')
#
# races_set = utils.RACES_SET
# count = 0
# for race in races_set:
#     my_df = pd.read_csv(f'data/real_data/races/{race}/id_amount')
#     count += len(my_df)
# print(f'check: {count}')

# my_dict = {'1': {'i': 1, 'j': 3},
#            '2': {'i': 5, 'j': 7}}
# print(list(my_dict.values()))

level = 'A'
race = 'AFA'
# df = pd.read_csv(f'data/real_data/levels/{level}/races/{race}/i_j_probability')
# sum_prob = 0.0
# for index, row in df.iterrows():
#     probability = row['probability']
#     sum_prob += probability
#
# print(sum_prob)

# df = pd.read_csv(f'data/real_data/levels/{level}/races/{race}/i_probability')
# sum_prob = 0.0
# for index, row in df.iterrows():
#     probability = row['probability']
#     sum_prob += probability
# print(sum_prob)

# total_probs = [
# 1.6359196801531626e-16,
# 4.63303574581979e-17,
# 6.818216970506675e-18,
# 3.733232977279645e-18,
# 3.0346568901841154e-18,
# 2.593806809530494e-18,
# 7.586644504163209e-19,
# 7.586644504163209e-19,
# 1.4740259796471218e-19
# ]
# total_sum = sum(total_probs)
#
# partial_prob214 = [
# 1.6359196801531626e-16,
# 3.733232977279645e-18
# ]
# partial_sum214 = sum(partial_prob214)
#
# partial_prob201 = [
# 4.63303574581979e-17,
# 6.818216970506675e-18,
# 3.0346568901841154e-18,
# 7.586644504163209e-19,
# 1.4740259796471218e-19
# ]
# partial_sum201 = sum(partial_prob201)
#
# partial_prob205 = [
# 2.593806809530494e-18,
# 7.586644504163209e-19
# ]
# partial_sum205 = sum(partial_prob205)
#
# print(f'total sum: {total_sum}')
# print(f' total sum for 214: {partial_sum214}')
# print(f' total sum for 201: {partial_sum201}')
# print(f' total sum for 205: {partial_sum205}')
#
# print(f' fraction for 214: {partial_sum214 / total_sum}')
# print(f' fraction for 201: {partial_sum201 / total_sum}')
# print(f' fraction for 205: {partial_sum205 / total_sum}')
# print(f'{partial_sum214 / total_sum + partial_sum201 / total_sum + partial_sum205 / total_sum}')

# plt.plot([0.0, 100.0], [0.0, 100.0])
#
# # scatter plot
# plt.scatter([3,4,5], [6,8,10], marker='.')
# plt.xlabel('Variance from formula')
# plt.ylabel('Variance over experiments')
# plt.show()

# z = np.array([[1,5,7], [3,8,9], [10,11,12]])
# couples = 3 * 4 // 2
# cdf = np.zeros(couples)
# index = 0
# for i in range(z.shape[0]):
#     for j in range(i, z.shape[1]):
#         cdf[index] = z[i, j]
#         if index > 0:
#             cdf[index] += cdf[index - 1]
#         index += 1
# print(cdf)
# print(len(cdf))

# alleles_count = 2
# population_amount = 3
# alpha_val = 1.0
# uncertainty_val = 0.0
# alleles_probabilities = np.array([0.5, 0.5])
#
# probabilities = np.zeros(shape=(alleles_count, alleles_count))
# for t in range(alleles_count):
#     for m in range(t, alleles_count):
#         if t == m:
#             probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (alleles_probabilities[m] ** 2)
#         else:
#             # we don't multiply by 2 here yet
#             probabilities[t, m] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
#             probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
# # matrix where every row is the alleles for a person
# alleles_individuals = np.zeros(
#     (population_amount, 2), dtype=np.int32)
#
# print(f'probabilities matrix:')
# print(probabilities)
#
# # 1) assignment of alleles with certainty
# for k in range(population_amount):
#     probabilities_list = probabilities.flatten()
#     print(f'probabilities flattened: {probabilities_list}')
#     index = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=1)[0]
#     print(f'index chosen: {index}')
#     col = index % alleles_count
#     row = (index - col) // alleles_count
#     print(f'alleles chosen: row: {row}, col: {col}')
#     row, col = min(row, col), max(row, col)
#     # assignment of the alleles to the person
#     alleles_individuals[k, :] = [row, col]
#
# print('assigned alleles:')
# print(alleles_individuals)
#
# # now we multiply the upper triangle by 2
# for t in range(alleles_count):
#     for m in range(t + 1, alleles_count):
#         probabilities[t, m] *= 2
#
# print(f'updated probabilities:')
# print(probabilities)
#
# # matrix p_k(i,j) = A[i,j,k]
# all_probabilities = np.zeros(shape=(alleles_count, alleles_count, population_amount))
#
# # adding uncertainty to our model
# for k in range(population_amount):
#     # person k has alleles j,l
#     j, l = alleles_individuals[k]
#     j, l = min(j, l), max(j, l)
#
#     # choice whether this person will have uncertain alleles
#     choice = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=1)[0]
#
#     # this person has certain alleles
#     if choice == 1:
#         all_probabilities[j, l, k] = 1.0
#     # this person has uncertain alleles
#     if choice == 0:
#         for t in range(alleles_count):
#             for m in range(t, alleles_count):
#                 all_probabilities[t, m, k] = probabilities[t, m]
#
# for i in range(population_amount):
#     print(f'probabilities of person {i}:')
#     print(all_probabilities[:,:,i])
#
#
# for t in range(alleles_count):
#     for m in range(t, alleles_count):
#         probability = 0.0
#         for k in range(population_amount):
#             probability += all_probabilities[t, m, k]
#         probability /= population_amount
#
#         probabilities[t, m] = probability
#
# observations = probabilities * population_amount
# print('observations:')
# print(observations)