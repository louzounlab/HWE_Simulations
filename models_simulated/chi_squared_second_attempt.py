import random
import numpy as np
import scipy.stats as stats
from matplotlib import pyplot as plt
import os
import csv


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
    return probs


# need to calculate only on the upper triangle because the matrices are symmetric
def calculate_chi_squared_value(alleles_amount, population_amount_, alleles_probabilities,
                                observed_probabilities, correction, cutoff,
                                test_type):
    value = 0.0
    amount_of_small_expected_ = 0
    for row in range(alleles_amount):
        for col in range(row, alleles_amount):
            mult = 1.0
            if row != col:
                mult = 2.0
            expected_val = mult * population_amount_ * alleles_probabilities[row] * alleles_probabilities[col]
            observed_val = population_amount_ * observed_probabilities[row, col]
            correction_val = correction[row, col]
            if test_type == 1:
                variance_val = expected_val * correction_val
            else:
                variance_val = expected_val
            if variance_val < cutoff:
                amount_of_small_expected_ += 1
                continue
            value += ((expected_val - observed_val) ** 2) / variance_val
    return value, amount_of_small_expected_


def generate_alleles_probabilities(alleles_count):
    return calculate_alleles_probabilities(alleles_count)


# generate simulated data as a list of lists (every row is: [k, i, j, p_k(i,j)])
def generate_data(alleles_count, population_amount, alleles_probabilities, alpha_val, uncertainty_val):
    # print('generate probablities')
    # probabilities {p(i)}
    # alleles_probabilities = calculate_alleles_probabilities(alleles_count)

    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (
                        alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
                probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]

    # matrix where every row is the alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, 2), dtype=np.int32)

    # print(f'alleles: {alleles_probabilities}')
    # print(f' marginal: {marginal_probabilities}')

    # 1) assignment of alleles with certainty
    probabilities_list = probabilities.flatten()
    # print('sample and append certain indices')
    indices = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=population_amount)
    for k in range(population_amount):
        # notice the alleles i != j have probability 2 * p(i,j)
        index = indices[k]
        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

    data = []
    # print('sample uncertain alleles')
    choices = random.choices(population=[0, 1], weights=[1 - uncertainty_val, uncertainty_val], k=population_amount)
    sum_uncertain_donors = sum(choices)
    observations_amount_for_uncertain_donor = 10
    indices_uncertain = random.choices(population=range(len(probabilities_list)),
                                       weights=probabilities_list,
                                       k=observations_amount_for_uncertain_donor * sum_uncertain_donors)
    current_uncertainty_index = 0
    # print('append uncertain alleles')
    for k in range(population_amount):
        # person k has alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)
        # no uncertainty
        if choices[k] == 0:
            data.append([k, j, l, 1.0])
        else:
            # matrix where each row represents alleles (i, j)
            alleles_uncertainty = np.zeros(shape=(observations_amount_for_uncertain_donor, 2), dtype=int)
            probabilities_uncertainty = np.zeros(shape=observations_amount_for_uncertain_donor)
            for t in range(observations_amount_for_uncertain_donor):
                index = indices_uncertain[current_uncertainty_index + t]
                # the right allele of this index element
                col = index % alleles_count
                # the left allele
                row = (index - col) // alleles_count
                # making sure i <= j
                row, col = min(row, col), max(row, col)
                # adding uncertain alleles to the person
                alleles_uncertainty[t, :] = [row, col]
                # getting the probability
                probability = probabilities[row, col]
                if row != col:
                    probability *= 2
                # adding the weight of the uncertain alleles
                probabilities_uncertainty[t] = probability

            # normalizing the uncertain observations of current person into probabilities
            probabilities_uncertainty = probabilities_uncertainty / np.sum(probabilities_uncertainty)

            # appending uncertain observations to current person
            for t in range(observations_amount_for_uncertain_donor):
                data.append([k,
                             alleles_uncertainty[t, 0],
                             alleles_uncertainty[t, 1],
                             probabilities_uncertainty[t]])
            current_uncertainty_index += observations_amount_for_uncertain_donor
    return data


def generate_data_old_simulation(alleles_count, population_amount, alpha_val, uncertainty_val,
                                 alleles_probabilities):
    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (
                        alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
                probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]

    # matrix where every row is the alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, 2), dtype=np.int32)

    # print(f'alleles: {alleles_probabilities}')
    # print(f' marginal: {marginal_probabilities}')

    # 1) assignment of alleles with certainty
    probabilities_list = probabilities.flatten()
    indices = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=population_amount)
    for k in range(population_amount):
        # notice the alleles i != j have probability 2 * p(i,j)
        index = indices[k]
        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

    # now we multiply the upper triangle by 2
    # for t in range(alleles_count):
    #     for m in range(t + 1, alleles_count):
    #         probabilities[t, m] *= 2

    # print(f'alleles_individuals: {alleles_individuals}')

    # matrix p_k(i,j) = A[i,j,k]
    # all_probabilities = np.zeros(shape=(alleles_count, alleles_count, population_amount))

    choices = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=population_amount)

    data = []
    # adding uncertainty to our model
    for k in range(population_amount):
        # person k has alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)

        # choice whether this person will have uncertain alleles
        choice = choices[k]

        # this person has certain alleles
        if choice == 1:
            # all_probabilities[j, l, k] = 1.0
            data.append([k, j, l, 1.0])
        # this person has uncertain alleles
        if choice == 0:
            # sum_ = 0
            for t in range(alleles_count):
                for u in range(t, alleles_count):
                    probability = probabilities[t, u]
                    if t != u:
                        probability *= 2
                    data.append([k, t, u, probability])
            # all_probabilities[:, :, k] = probabilities
    return data


# given data: list of lists, run Traditional Chi Squared and Asta and return both significance result.
def run_experiment(data, cutoff_value=0.0):
    id_to_index = {}
    allele_to_index = {}
    index_to_allele = {}
    # index1_index2 -> [obs_prob, corr]

    # id -> {'i_j' -> [i, j, O_k(i,j)], ...}
    # actual id string, i and j are indices
    id_to_i_j_to_i_j_observation = {}

    # first read all the rows and get indices of ids and alleles and amounts
    for lst in data:
        id_row = str(lst[0])
        allele_1 = str(lst[1])
        allele_2 = str(lst[2])
        # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
        observation = float(lst[3])

        if id_row not in id_to_index:
            id_to_index[id_row] = len(id_to_index)

        if allele_1 not in allele_to_index:
            # updating the inverse dictionary
            index_to_allele[len(allele_to_index)] = allele_1
            # updating the dictionary
            allele_to_index[allele_1] = len(allele_to_index)
        if allele_2 not in allele_to_index:
            # updating the inverse dictionary
            index_to_allele[len(allele_to_index)] = allele_2
            # updating the dictionary
            allele_to_index[allele_2] = len(allele_to_index)

        # get indices of alleles and make i<=j
        i = allele_to_index[allele_1]
        j = allele_to_index[allele_2]
        i, j = min(i, j), max(i, j)
        i_j = f'{i}_{j}'
        # update dictionary
        if id_row not in id_to_i_j_to_i_j_observation:
            id_to_i_j_to_i_j_observation[id_row] = {i_j: [i, j, observation]}
        else:
            # id exists in dictionary, we need to append the observation for i, j
            # if i_j exists, we need to add the probability and if it doesn't,
            # we need to append a new observation

            # if observation i, j doesn't exist
            if i_j not in id_to_i_j_to_i_j_observation[id_row]:
                id_to_i_j_to_i_j_observation[id_row][i_j] = [i, j, observation]
            else:
                # observation i, j exists. add the new observation
                id_to_i_j_to_i_j_observation[id_row][i_j][2] += observation

    alleles_count = len(allele_to_index)
    population_amount = len(id_to_index)

    # {p(i)}
    alleles_probabilities = np.zeros(alleles_count)

    # # {p(i, j)}
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))
    #
    # # correction matrix
    correction = np.zeros(shape=(alleles_count, alleles_count))

    # go over the ids
    for current_id in id_to_i_j_to_i_j_observation:
        # sum observations
        sum_observations = 0.0
        for i_j in id_to_i_j_to_i_j_observation[current_id]:
            sum_observations += id_to_i_j_to_i_j_observation[current_id][i_j][2]
        # get the probabilities
        for i_j in id_to_i_j_to_i_j_observation[current_id]:
            # i and j are already i <= j
            i = id_to_i_j_to_i_j_observation[current_id][i_j][0]
            j = id_to_i_j_to_i_j_observation[current_id][i_j][1]
            probability = id_to_i_j_to_i_j_observation[current_id][i_j][2] / sum_observations

            alleles_probabilities[i] += 0.5 * probability
            alleles_probabilities[j] += 0.5 * probability
            observed_probabilities[i, j] += probability
            correction[i, j] += (probability ** 2)

    # p(i) = sum_k_j p_k(i,j) / N
    alleles_probabilities /= population_amount

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            observed_probabilities[i, j] /= population_amount
            if observed_probabilities[i, j] == 0:
                correction[i, j] = 1.0
            else:
                correction[i, j] /= (population_amount * observed_probabilities[i, j])

    chi_squared_stat_old, amount_of_small_expected_old = calculate_chi_squared_value(alleles_amount=alleles_count,
                                                                                     population_amount_=population_amount,
                                                                                     alleles_probabilities=alleles_probabilities,
                                                                                     observed_probabilities=observed_probabilities,
                                                                                     correction=correction,
                                                                                     cutoff=cutoff_value,
                                                                                     test_type=0)
    couples_amount = (alleles_count * (alleles_count + 1)) / 2 - 1
    dof_old = couples_amount - amount_of_small_expected_old

    p_value_old = 1 - stats.chi2.cdf(x=chi_squared_stat_old,
                                     df=dof_old)

    chi_squared_stat_new, amount_of_small_expected_new = calculate_chi_squared_value(alleles_amount=alleles_count,
                                                                                     population_amount_=population_amount,
                                                                                     alleles_probabilities=alleles_probabilities,
                                                                                     observed_probabilities=observed_probabilities,
                                                                                     correction=correction,
                                                                                     cutoff=cutoff_value,
                                                                                     test_type=1)
    dof_new = couples_amount - amount_of_small_expected_new

    p_value_new = 1 - stats.chi2.cdf(x=chi_squared_stat_new,
                                     df=dof_new)
    # return int(p_value_old < 0.05), int(p_value_new < 0.05)
    return p_value_old, p_value_new, dof_old, dof_new


# given data: list of lists, run Traditional Chi Squared and Asta and return both significance result.
def plot_variance_vs_corrected_variance(data):
    id_to_index = {}
    allele_to_index = {}
    index_to_allele = {}
    # index1_index2 -> [obs_prob, corr]

    # id -> {'i_j' -> [i, j, O_k(i,j)], ...}
    # actual id string, i and j are indices
    id_to_i_j_to_i_j_observation = {}

    # first read all the rows and get indices of ids and alleles and amounts
    for lst in data:
        id_row = str(lst[0])
        allele_1 = str(lst[1])
        allele_2 = str(lst[2])
        # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
        observation = float(lst[3])

        if id_row not in id_to_index:
            id_to_index[id_row] = len(id_to_index)

        if allele_1 not in allele_to_index:
            # updating the inverse dictionary
            index_to_allele[len(allele_to_index)] = allele_1
            # updating the dictionary
            allele_to_index[allele_1] = len(allele_to_index)
        if allele_2 not in allele_to_index:
            # updating the inverse dictionary
            index_to_allele[len(allele_to_index)] = allele_2
            # updating the dictionary
            allele_to_index[allele_2] = len(allele_to_index)

        # get indices of alleles and make i<=j
        i = allele_to_index[allele_1]
        j = allele_to_index[allele_2]
        i, j = min(i, j), max(i, j)
        i_j = f'{i}_{j}'
        # update dictionary
        if id_row not in id_to_i_j_to_i_j_observation:
            id_to_i_j_to_i_j_observation[id_row] = {i_j: [i, j, observation]}
        else:
            # id exists in dictionary, we need to append the observation for i, j
            # if i_j exists, we need to add the probability and if it doesn't,
            # we need to append a new observation

            # if observation i, j doesn't exist
            if i_j not in id_to_i_j_to_i_j_observation[id_row]:
                id_to_i_j_to_i_j_observation[id_row][i_j] = [i, j, observation]
            else:
                # observation i, j exists. add the new observation
                id_to_i_j_to_i_j_observation[id_row][i_j][2] += observation

    alleles_count = len(allele_to_index)
    population_amount = len(id_to_index)

    # {p(i)}
    alleles_probabilities = np.zeros(alleles_count)

    # # {p(i, j)}
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))
    #
    # # correction matrix
    correction = np.zeros(shape=(alleles_count, alleles_count))

    # go over the ids
    for current_id in id_to_i_j_to_i_j_observation:
        # sum observations
        sum_observations = 0.0
        for i_j in id_to_i_j_to_i_j_observation[current_id]:
            sum_observations += id_to_i_j_to_i_j_observation[current_id][i_j][2]
        # get the probabilities
        for i_j in id_to_i_j_to_i_j_observation[current_id]:
            # i and j are already i <= j
            i = id_to_i_j_to_i_j_observation[current_id][i_j][0]
            j = id_to_i_j_to_i_j_observation[current_id][i_j][1]
            probability = id_to_i_j_to_i_j_observation[current_id][i_j][2] / sum_observations

            alleles_probabilities[i] += 0.5 * probability
            alleles_probabilities[j] += 0.5 * probability
            observed_probabilities[i, j] += probability
            correction[i, j] += (probability ** 2)

    # p(i) = sum_k_j p_k(i,j) / N
    alleles_probabilities /= population_amount

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            observed_probabilities[i, j] /= population_amount
            if observed_probabilities[i, j] == 0:
                correction[i, j] = 1.0
            else:
                correction[i, j] /= (population_amount * observed_probabilities[i, j])
    expected_vecter = []
    expected_corrected_vector = []
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            expected = population_amount * alleles_probabilities[i] * alleles_probabilities[j]
            if i != j:
                expected *= 2
            expected_corrected = expected * correction[i, j]
            if expected_corrected < 2.0:
                continue
            expected_vecter.append(expected)
            expected_corrected_vector.append(expected_corrected)
    plt.scatter(expected_corrected_vector, expected_vecter, color='deeppink')
    max_val = max(expected_corrected_vector + expected_vecter)
    min_val = min(expected_corrected_vector + expected_vecter)
    plt.xlim([min_val, max_val])
    plt.ylim([min_val, max_val])
    plt.xlabel('Expected * Corrected')
    plt.ylabel('Expected')
    plt.show()

