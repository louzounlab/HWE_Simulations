import random
import numpy as np
import scipy.stats as stats
import os
import csv
import time

# # assuming x is 2-dimensional
# def softmax_1d(x):
#     sum = 0.0
#     for row in range(x.shape[0]):
#         sum += np.exp(x[row])
#     return np.exp(x) / sum
#
#
# # generate a vector of probabilities
# def calculate_alleles_probabilities(alleles_count):
#     probs = np.random.uniform(0.0, 2.0, size=alleles_count)
#     probs = softmax_1d(probs)
#     return probs
#
#
# # returns the index of element closest to target.
# # cdf=[p1, p1+p2, ..., 1]
# def binary_search(cdf, target: float = 0):
#     # for i in range(len(cdf)):
#     #     if cdf[i][0] >= target:
#     #         return i
#     start = 0
#     end = len(cdf)
#     while start < end:
#         mid = (end + start) // 2
#         # print(f'start: {start}. mid: {mid}. end: {end}')
#         if cdf[mid] < target:
#             start = mid + 1
#         else:
#             end = mid
#     # now start and end pointing to the elements closest to 0
#     # pick the index of the closer one
#     return start
#
#
# def get_cdf(probabilities):
#     cdf = np.zeros(shape=len(probabilities))
#     current_sum = 0.0
#     for i in range(len(probabilities)):
#         current_sum += probabilities[i]
#         cdf[i] = current_sum
#     return cdf
#
#
# # given cdf and k (amount to sample), sample k indices from cdf
# def sample_from_cdf(cdf, k=1):
#     indices = []
#     uniforms = np.random.uniform(0, 1, size=k)  # sample k U[0,1) random variables
#     for i in range(k):
#         index = binary_search(cdf=cdf, target=uniforms[i])
#         indices.append(index)
#     return indices
#
#
# alleles_count = 1000
# population_amount = 1000000
# alpha_val = 0.96
# uncertainty_val = 0.2
# # probabilities {p(i)}
# # print('calculating {p(i)}')
# alleles_probabilities = calculate_alleles_probabilities(alleles_count)
# # print('calculating {p(i,j)}')
# probabilities = np.zeros(shape=(alleles_count, alleles_count))
# for t in range(alleles_count):
#     for m in range(t, alleles_count):
#         if t == m:
#             probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (
#                     alleles_probabilities[m] ** 2)
#         else:
#             # we don't multiply by 2 here yet
#             probabilities[t, m] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
#             probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
#
# # print('sampling alleles')
# # matrix where every row is the alleles for a person
# alleles_individuals = np.zeros(
#     (population_amount, 2), dtype=np.int32)
#
# # print(f'alleles: {alleles_probabilities}')
# # print(f' marginal: {marginal_probabilities}')
#
# # 1) assignment of alleles with certainty
# probabilities_list = probabilities.flatten()
# cdf_list = get_cdf(probabilities_list)
# print('sampling indices')
# start_time = time.time()
# # for each donor, sample (i, j) according to {p(i,j)}
# indices = sample_from_cdf(cdf=cdf_list, k=population_amount)
# end_time = time.time()
# print(f'time elapsed: {end_time - start_time}')
# print('sampling indices using random')
# start_time = time.time()
# choices = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=population_amount)
# end_time = time.time()
# print(f'time elapsed: {end_time - start_time}')
# print('appending incides')
# for k in range(population_amount):
#     # notice the alleles i != j have probability 2 * p(i,j)
#     index = indices[k]
#     # the right allele of this index element
#     col = index % alleles_count
#     # the left allele
#     row = (index - col) // alleles_count
#     # making sure i <= j
#     row, col = min(row, col), max(row, col)
#     # assignment of the alleles to the person
#     alleles_individuals[k, :] = [row, col]

# print('calculating uncertainties')
# choices = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=population_amount)
# print('adding uncertainties')
# data = []
# # adding uncertainty to our model
# for k in range(population_amount):
#     # person k has alleles j,l
#     j, l = alleles_individuals[k]
#     j, l = min(j, l), max(j, l)
#
#     # choice whether this person will have uncertain alleles
#     choice = choices[k]
#
#     # this person has certain alleles
#     if choice == 1:
#         # all_probabilities[j, l, k] = 1.0
#         # writer.writerow([k, j, l, 1.0])
#         data.append([k, j, l, 1.0])
#     # this person has uncertain alleles
#     if choice == 0:
#         # all_probabilities[:, :, k] = probabilities
#         # sample 10 indices from cdf
#         indices_uncertainty = sample_from_cdf(cdf=cdf_list, k=10)
#         # matrix where each row represents alleles (i, j)
#         alleles_uncertainty = np.zeros(shape=(len(indices_uncertainty), 2))
#         probabilities_uncertainty = np.zeros(shape=len(indices_uncertainty))
#         for i in range(len(indices_uncertainty)):
#             index = indices[i]
#             # the right allele of this index element
#             col = index % alleles_count
#             # the left allele
#             row = (index - col) // alleles_count
#             # making sure i <= j
#             row, col = min(row, col), max(row, col)
#             # adding uncertain alleles to the person
#             alleles_uncertainty[i, :] = [row, col]
#             # getting the probability
#             probability = probabilities[row, col]
#             if row != col:
#                 probability *= 2
#             # adding the weight of the uncertain alleles
#             probabilities_uncertainty[i] = probability
#
#         # normalizing the uncertain observations of current person into probabilities
#         probabilities_uncertainty = probabilities_uncertainty / np.sum(probabilities_uncertainty)
#
#         # appending uncertain observations to current person
#         for i in range(len(indices_uncertainty)):
#             data.append([k,
#                          alleles_uncertainty[i, 0],
#                          alleles_uncertainty[i, 1],
#                          probabilities_uncertainty[i]])


def index_to_indices(idx, mat):
    rows = mat.shape[0]
    cols = mat.shape[1]

    col = idx % cols
    row = (idx - col) // rows
    print(f'{row}, {col}')


z = np.array([[4, 5, 6], [1, 8, 3]])
z_flat = z.flatten()
print(len(z_flat))
print(z_flat)
print((-z_flat).argsort())
print(f'largest index: {(-z_flat).argsort()[0]}')
print('largest indices:')

index_to_indices((-z_flat).argsort()[0], z)
