import random
import numpy as np
import scipy.stats as stats
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


# returns the index of element closest to target.
# cdf=[p1, p1+p2, ..., 1]
def binary_search(cdf, target: float = 0):
    # for i in range(len(cdf)):
    #     if cdf[i][0] >= target:
    #         return i
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


def get_cdf(probabilities):
    cdf = np.zeros(shape=len(probabilities))
    current_sum = 0.0
    for i in range(len(probabilities)):
        current_sum += probabilities[i]
        cdf[i] = current_sum
    return cdf


# given cdf and k (amount to sample), sample k indices from cdf
def sample_from_cdf(cdf, k=1):
    indices = []
    uniforms = np.random.uniform(0, 1, size=k)  # sample k U[0,1) random variables
    for i in range(k):
        index = binary_search(cdf=cdf, target=uniforms[i])
        indices.append(index)
    return indices


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


# generate simulated data as a list of lists (every row is: [k, i, j, p_k(i,j)])
def generate_data(alleles_count, population_amount, alpha_val, uncertainty_val):
    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)

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
    cdf_list = get_cdf(probabilities_list)

    # for each donor, sample (i, j) according to {p(i,j)}
    indices = sample_from_cdf(cdf=cdf_list, k=population_amount)
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

    # if os.path.exists(file_name) and os.path.isfile(file_name):
    #     os.remove(file_name)
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
            # writer.writerow([k, j, l, 1.0])
            data.append([k, j, l, 1.0])
        # this person has uncertain alleles
        if choice == 0:
            # all_probabilities[:, :, k] = probabilities
            # sample 10 indices from cdf
            indices_uncertainty = sample_from_cdf(cdf=cdf_list, k=10)
            # matrix where each row represents alleles (i, j)
            alleles_uncertainty = np.zeros(shape=(len(indices_uncertainty), 2))
            probabilities_uncertainty = np.zeros(shape=len(indices_uncertainty))
            for i in range(len(indices_uncertainty)):
                index = indices[i]
                # the right allele of this index element
                col = index % alleles_count
                # the left allele
                row = (index - col) // alleles_count
                # making sure i <= j
                row, col = min(row, col), max(row, col)
                # adding uncertain alleles to the person
                alleles_uncertainty[i, :] = [row, col]
                # getting the probability
                probability = probabilities[row, col]
                if row != col:
                    probability *= 2
                # adding the weight of the uncertain alleles
                probabilities_uncertainty[i] = probability

            # normalizing the uncertain observations of current person into probabilities
            probabilities_uncertainty = probabilities_uncertainty / np.sum(probabilities_uncertainty)

            # appending uncertain observations to current person
            for i in range(len(indices_uncertainty)):
                data.append([k,
                             alleles_uncertainty[i, 0],
                             alleles_uncertainty[i, 1],
                             probabilities_uncertainty[i]])
            #
            # for y in range(alleles_count):
            #     for u in range(y, alleles_count):
            #         if y == u:
            #             writer.writerow([k, y, u, probabilities[y, u]])
            #         else:
            #             # probabilities matrix is symmetric, not upper triangular
            #             writer.writerow([k, y, u, 2 * probabilities[y, u]])
    return data


# given data: list of lists, run Traditional Chi Squared and Asta and return both significance result.
def run_experiment(data, cutoff_value=0.0):
    id_to_index = {}
    allele_to_index = {}
    index_to_allele = {}
    # index1_index2 -> [obs_prob, corr]

    # first read all the rows and get indices of ids and alleles and amounts
    for lst in data:
        current_id = lst[0]
        allele_1 = lst[1]
        allele_2 = lst[2]
        # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
        # probability = float(lst[3])

        if current_id not in id_to_index:
            id_to_index[current_id] = len(id_to_index)

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

    alleles_count = len(allele_to_index)
    population_amount = len(id_to_index)

    # {p(i)}
    alleles_probabilities = np.zeros(alleles_count)

    # # {p(i, j)}
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))
    #
    # # correction matrix
    correction = np.zeros(shape=(alleles_count, alleles_count))

    # calculate {p_k(i,j)}
    for lst in data:
        # id = lst[0]
        allele_1 = lst[1]
        allele_2 = lst[2]
        # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
        probability = float(lst[3])

        # id_index = id_to_index[id]

        allele_1_index = allele_to_index[allele_1]
        allele_2_index = allele_to_index[allele_2]

        allele_1_index, allele_2_index = min(allele_1_index, allele_2_index), max(allele_1_index, allele_2_index)

        alleles_probabilities[allele_1_index] += 0.5 * probability
        alleles_probabilities[allele_2_index] += 0.5 * probability

        observed_probabilities[allele_1_index, allele_2_index] += probability

        correction[allele_1_index, allele_2_index] += (probability ** 2)

    # p(i) = sum_k_j p_k(i,j) / N
    alleles_probabilities /= population_amount

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            observed_probabilities[i, j] /= population_amount
            if observed_probabilities[i, j] == 0:
                correction[i, j] = 1.0
            else:
                correction[i, j] /= (population_amount * observed_probabilities[i, j])

    # for i in range(alleles_count):
    #     for j in range(alleles_count):
    #         observed_probabilities[i, j] /= population_amount
    #         if observed_probabilities[i, j] == 0:
    #             correction[i, j] = 1.0
    #         else:
    #             correction[i, j] /= (population_amount * observed_probabilities[i, j])

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

    return int(p_value_old < 0.05), int(p_value_new < 0.05)
