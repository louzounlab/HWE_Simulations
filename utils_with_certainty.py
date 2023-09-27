import numpy as np
import random
import array as arr
import math
import matplotlib.pyplot as plt


###################
# UTILITY FUNCTIONS
###################


# assuming x is 2-dimensional
def softmax_1d(x):
    sum = 0.0
    for row in range(x.shape[0]):
        sum += np.exp(x[row])
    return np.exp(x) / sum


# generate a vector of probabilities
def calculate_alleles_probabilities(alleles_count):
    probs = np.random.uniform(0.0, 2.0, size=alleles_count)
    probs = softmax_1d(probs)
    # probs = np.zeros(alleles_count) + (1 / alleles_count)
    return probs


# generate a marginal probabilities matrix
# column i represent the probabilities p_i(j), j are for the rows
def calculate_marginal_probabilities(alleles_count):
    probs = np.zeros(shape=(alleles_count, alleles_count))
    for i_ in range(alleles_count):
        marginal = calculate_alleles_probabilities(alleles_count)
        probs[:, i_] = marginal
    return probs


# data: N x 2, row = (j, k). alleles j,k
def count_row_occurrence_in_2d_array(j, k, data):
    count = 0
    # for every row (person) in population
    for row in range(data.shape[0]):
        if (data[row, 0] == j and data[row, 1] == k) or (data[row, 0] == k and data[row, 1] == j):
            count += 1
            # count += (alleles_probabilities[j] * alleles_probabilities[k] / (population_probs[j,k,row]))
            # count += 1

    return count


def calculate_probability_given_i_j(marginal_probabilities_, i_, j_, k_1, k_2):
    prob = marginal_probabilities_[k_1, i_] * marginal_probabilities_[k_2, j_]
    # EVEN
    if i_ != j_:
        prob += marginal_probabilities_[k_1, j_] * marginal_probabilities_[k_2, i_]
    return prob


def calc_O_ij(population_amount, allele_count, alpha_val, alleles_probabilities, probabilities):
    observed = np.zeros(shape=(allele_count, allele_count))
    for i in range(allele_count):
        for j in range(i, allele_count):
            if i == j:
                observed[i, j] = np.random.binomial(population_amount, get_p_ij(probabilities, i, j) * (1 + alpha_val))
            else:
                observed[i, j] = np.random.binomial(
                    population_amount, get_p_ij(probabilities, i, j))
    return observed


def calc_O_ij_test(population_amount, allele_count, alpha_val, alleles_probabilities):
    observed = np.zeros(shape=(allele_count, allele_count))
    for i in range(allele_count):
        for j in range(i, allele_count):
            if i == j:
                observed[i, j] = 1000
            else:
                observed[i, j] = 0
    return observed


# calc matrix of {p(i,j)}
def calc_p_ij(alleles_count, alleles_probabilities):
    probs = np.zeros(shape=(alleles_count, alleles_count))

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            if i == j:
                probs[i, j] = (alleles_probabilities[i] ** 2)
            else:
                probs[i, j] = 2 * alleles_probabilities[i] * alleles_probabilities[j]

    return probs


# instead of using a {p(i,j)} matrix, we can use this function to access the elements
def get_p_ij(probabilities, i, j):
    i, j = min(i, j), max(i, j)
    return probabilities[i, j]


# returns O(i,j) from the upper triangle. order of i,j doesnt matter
def get_O_ij(observed, i, j):
    row, col = min(i, j), max(i, j)
    return observed[row, col]


# calculate a dictionary from alleles to cdfs: (cdf_dict is the dictionary)
# allele i -> [(O(0,i)/n, 0), (O(0,i) + O(1,i) / n,1), ...,(1,n)]
# every value for key i is a list of tuples, left element is the cdf of observed row,
# right element in the tuple is the relevant index that goes up
def calculate_cdf_dict(alleles_count, observed, cdf_dict):
    # cdf_dict is the dictionary in which we will calculate the cdf,
    # we need to remove all the items before.
    cdf_dict.clear()
    # go over the keys
    for i in range(alleles_count):
        # this list will be the value for the key i
        cdf_list = []
        # first calculate the sum of {O(i,k)}_k
        sum_row = 0
        # we multiply by half outside the diagonal, because without it, the probability is multiplied by 2
        for k in range(alleles_count):
            mult = 0.5
            if k == i:
                mult = 1
            sum_row += mult * get_O_ij(observed, i, k)
        # now fill the list
        current_sum = 0
        # now we go over the alleles, first when k<=i we will take min(i,k) and after that max(i,k)
        for k in range(alleles_count):
            mult = 0.5
            if k == i:
                mult = 1
            current_sum += mult * get_O_ij(observed, i, k)
            # # given allele i, we dont want to choose it as the second allele
            # if k != i:
            cdf_list.append((current_sum / sum_row, k))
        cdf_dict[i] = cdf_list


# given a row - i, we will update the row based on the observations
def update_row_in_cdf_dict(alleles_count, observed, cdf_dict, i):
    cdf_list = []
    # first calculate the sum of {O(i,k)}_k
    sum_row = 0
    # we multiply by half outside the diagonal, because without it, the probability is multiplied by 2
    for k in range(alleles_count):
        mult = 0.5
        if k == i:
            mult = 1
        sum_row += mult * get_O_ij(observed, i, k)
    # now fill the list
    current_sum = 0
    # now we go over the alleles, first when k<=i we will take min(i,k) and after that max(i,k)
    for k in range(alleles_count):
        mult = 0.5
        if k == i:
            mult = 1
        current_sum += mult * get_O_ij(observed, i, k)
        # # given allele i, we dont want to choose it as the second allele
        # if k != i:
        cdf_list.append((current_sum / sum_row, k))
    cdf_dict[i] = cdf_list


# returns [(p_00, 0), (p_00+p_01, 1), ..., (p_00+...+p_kk, k)]
def calculate_cdf(alleles_count, probabilities):
    size = (alleles_count * (alleles_count + 1)) // 2
    cdf_probabilities = arr.array('f', [2 for i in range(size)])
    cdf_alleles = arr.array('i', [2 for i in range(size)])
    # lst = []
    current_sum = 0
    index = 0
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            current_sum += probabilities[i, j]
            # lst.append((current_sum, j))
            cdf_probabilities[index] = current_sum
            cdf_alleles[index] = j
            index += 1
    return cdf_probabilities, cdf_alleles


# returns the index of element closest to target.
# cdf=[(p_11, 1), (p_11+p_12, 2), ..., (p_11+...+p_kk, k)]
def binary_search(cdf, target: float = 0):
    # for i in range(len(cdf)):
    #     if cdf[i][0] >= target:
    #         return i
    start = 0
    end = len(cdf)
    while start < end:
        mid = (end + start) // 2
        # print(f'start: {start}. mid: {mid}. end: {end}')
        if cdf[mid][0] < target:
            start = mid + 1
        else:
            end = mid
    # now start and end pointing to the elements closest to 0
    # pick the index of the closer one
    return start


def get_allele_from_cdf_dict(alleles_count, cdf_dict, observed, left_alleles: list):
    left_allele = left_alleles[1]
    # update row in dict
    update_row_in_cdf_dict(alleles_count, observed, cdf_dict, left_allele)

    random_probability = np.random.uniform(0, 1)
    index = binary_search(cdf_dict[left_allele], random_probability)
    # accessing the dict, then the list and then the right element of tuplef
    right_allele = cdf_dict[left_allele][index][1]

    # counter = 0
    counter_for_recalculating_cdf_dict = 0

    if get_O_ij(observed, left_alleles[1], right_allele) <= 0:
        raise ValueError('A very specific bad thing happened.')


        # o_ij = observed[min(left_allele, right_allele), max(left_allele, right_allele)]
        # print(f'O({left_allele},{right_allele}) = {o_ij}')
        # print(f'left allele: {left_alleles[1]}. right allele: {right_allele}. possibles O({left_allele},?):')
        # debug_list = []
        # for k in range(len(cdf_dict[left_allele])):
        #     debug_list.append(get_O_ij(observed, left_alleles[1], k))
        # print(debug_list)
        # print('dict for left allele:')
        # print(cdf_dict[left_allele])


        # left_allele = left_alleles[1]
        # random_probability = np.random.uniform(0, 1)
        # index = binary_search(cdf_dict[left_allele], random_probability)
        # # accessing the dict, then the list and then the right element of tuple
        # right_allele = cdf_dict[left_allele][index][1]

        # counter += 1
        # counter_for_recalculating_cdf_dict += 1
        # if counter >= 20:
        #     left_alleles[0], left_alleles[1] = left_alleles[1], left_alleles[0]
        #     counter = 0

        # if this happens, the cdf dictionary is too old, and we need to calculate it again.
        # if counter_for_recalculating_cdf_dict >= 20000:
        #     print('INITIALIZING CDF DICTIONARY!!!!!!!!!!!!!!!!!!!!!')
        #     calculate_cdf_dict(alleles_count, observed, cdf_dict=cdf_dict)
        #     counter_for_recalculating_cdf_dict = 0

    return right_allele


# cdf=[(p_11, 1), (p_11+p_12, 2), ..., (p_11+...+p_kk, k)]
# observed O(i,j)
# random_probability between 0 and 1
# left_allele: i, we will pick an allele j from the cdf such that O(i,j) > 0.
# returns the k corresponding to the closest value p_11+...+p_kk to the given probability
def get_allele_from_cdf(cdf_probabilities, cdf_alleles, observed, left_alleles: list):
    # debug_val = 1
    # # debug if there is no possibility to choose right allele such that the observed > 1:
    # for allele in range(observed.shape[0]):
    #     x = get_O_ij(observed, left_allele, allele)
    #     if x > 0:
    #         debug_val = 0
    # if debug_val == 1:
    #     print(f'cant choose an allele')
    # else:
    #     print('can')

    random_probability = np.random.uniform(0, 1)
    # convert cdf list to a list of probabilities
    # arr = [cdf[i][0] for i in range(len(cdf))]
    # take the index with closest value to random_probability
    index = binary_search(cdf_probabilities, random_probability)
    right_allele = cdf_alleles[index]

    counter = 0

    while get_O_ij(observed, left_alleles[0], right_allele) <= 0:
        random_probability = np.random.uniform(0, 1)
        # convert cdf list to a list of probabilities
        # take the index with closest value to random_probability
        index = binary_search(cdf_probabilities, random_probability)
        right_allele = cdf_alleles[index]
        # print('now')
        counter += 1
        if counter >= 20:
            left_alleles[0], left_alleles[1] = left_alleles[1], left_alleles[0]
            counter = 0
        #     sum_observed = 0
        #     for i in range(observed.shape[0]):
        #         for j in range(i, observed.shape[1]):
        #             sum_observed += observed[i, j]
        #     print(f'sum observed: {sum_observed}')
    return right_allele


# given O_ij, we pick two alleles (i,j) such that O(i,j) > 0
def calculate_couple_of_alleles(alleles_count, observed):
    couple = random.choices(population=range(alleles_count), k=2)
    couple[0], couple[1] = min(couple[0], couple[1]), max(couple[0], couple[1])
    # couple[0], couple[1] = int(couple[0]), int(couple[1])
    counter = 0
    while observed[couple[0], couple[1]] <= 0:
        couple = random.choices(population=range(alleles_count), k=2)
        couple[0], couple[1] = min(couple[0], couple[1]), max(couple[0], couple[1])
        # couple[0], couple[1] = int(couple[0]), int(couple[1])

        # we might have an infinite loop
        counter += 1
        if counter >= 100:
            print('All O_ij are probably zeros !!!')
            counter = 0

    return couple


# takes two indices i, j (order doesnt matter) and perform: O(i,j) += value_to_add.
# we only care for the upper triangle anyway.
# Python takes the reference of the array so the change is inplace.
# couple_indices = [i, j]
def modify_observed_ij(observed, couple_indices, value_to_add):
    i, j = min(couple_indices[0], couple_indices[1]), max(couple_indices[0], couple_indices[1])
    observed[i, j] += value_to_add
    # print(f'modify O({i}, {j}) + {value_to_add}')


# here we calculate ln(p_i) (the delta_i). given the observed and probabilities.
# couples is a matrix k x 3 where k is the amount of couples we take in the calculation.
# and every row is the couple i,j and a value of 1 or -1.
# population_amount_calculated represents the calculated amount of population based on observations.
# def update_current_delta_probability(list_probabilities: list, observed, probabilities, couples,
#                                      population_amount_calculated):
#     sum_current = 0
#     for row in range(couples.shape[0]):
#         # here we make sure that i <= j so we can access the upper triangle of probabilities
#         i, j = min(couples[row, 0], couples[row, 1]), max(couples[row, 0], couples[row, 1])
#         # i, j = int(i), int(j)
#
#         # z = np.log(population_amount_calculated * probabilities[i, j] / (1 - probabilities[i, j]))
#         sign = couples[row, 2]
#         val = 0
#         if sign == -1:
#             val = np.log(observed[i, j]) - np.log(population_amount_calculated - observed[i, j] + 1) \
#                   - np.log(probabilities[i, j]) + np.log(1 - probabilities[i, j])
#         elif sign == 1:
#             val = np.log(population_amount_calculated - observed[i, j]) - np.log(observed[i, j]) \
#                   + np.log(probabilities[i, j]) - np.log(1 - probabilities[i, j])
#
#         # print(observed)
#         # print(f'i: {i}, j: {j}')
#         # sum_current += sign * (z + np.log(observed[i, j]))
#         sum_current += val
#     list_probabilities.append(sum_current)
#
#     for row in range(couples.shape[0]):
#         # here we make sure that i <= j so we can access the upper triangle of probabilities
#         i, j = min(couples[row, 0], couples[row, 1]), max(couples[row, 0], couples[row, 1])
#         modify_observed_ij(observed=observed, couple_indices=[i, j],
#                            value_to_add=-1)


def update_current_delta_probability(list_probabilities: list, observed, alleles_probabilities, couples,
                                     population_amount_calculated, probabilities):
    for row in range(couples.shape[0]):
        # here we make sure that i <= j so we can access the upper triangle of probabilities
        i, j = min(couples[row, 0], couples[row, 1]), max(couples[row, 0], couples[row, 1])
        # i, j = int(i), int(j)

        # z = np.log(population_amount_calculated * probabilities[i, j] / (1 - probabilities[i, j]))
        sign = couples[row, 2]
        val = 0
        if sign == -1:
            # val = np.log(observed[i, j] / (population_amount_calculated - observed[i, j] + 1)) \
            #       + np.log((1 - get_p_ij(alleles_probabilities, i, j)) / get_p_ij(alleles_probabilities, i, j))
            val = np.log(get_O_ij(observed, i, j)) - np.log(population_amount_calculated - get_O_ij(observed, i, j) + 1) \
                  - np.log(get_p_ij(probabilities, i, j)) + np.log(1 - get_p_ij(probabilities, i, j))
        elif sign == 1:
            val = np.log(population_amount_calculated - get_O_ij(observed, i, j)) - np.log(get_O_ij(observed, i, j) + 1) \
                  + np.log(get_p_ij(probabilities, i, j)) - np.log(1 - get_p_ij(probabilities, i, j))


        # print(observed)
        # print(f'i: {i}, j: {j}')
        # sum_current += sign * (z + np.log(observed[i, j]))
        list_probabilities.append(val)
        # unbalanced_alleles[i] += sign
        # unbalanced_alleles[j] += sign
        modify_observed_ij(observed=observed, couple_indices=[i, j],
                           value_to_add=sign)


# calculate the starting index for the gibbs sampling
def calculate_start_time(alleles_count, population_amount, alleles_probabilities, observed, cdf_dict):
    start_time = 0.0
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            mult = 1
            if i != j:
                mult = 2
            val = observed[i, j] - mult * population_amount * alleles_probabilities[i] * alleles_probabilities[j]
            if val >= 0:
                val = abs(val) * cdf_dict[i][j][0] * (1 - alleles_probabilities[j])
            else:
                val = abs(val) * (1 - cdf_dict[i][j][0]) * alleles_probabilities[j]
            start_time += val
    return int(start_time)

    pass


def debug_observed(observed):
    non_zero_counter = 0
    for i in range(observed.shape[0]):
        for j in range(i, observed.shape[1]):
            if observed[i, j] <= 0:
                non_zero_counter += 1
    print(f'non zero elements in upper triangle of observed: {non_zero_counter}')
