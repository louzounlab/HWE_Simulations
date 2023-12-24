import random
import numpy as np
import scipy.stats as stats
from matplotlib import pyplot as plt


def softmax_1d(x):
    """
    Performs a softmax on a given numpy array x
    :param x: 1 dimensional numpy array
    :return: probabilities
    """
    sum = 0.0
    for row in range(x.shape[0]):
        sum += np.exp(x[row])
    return np.exp(x) / sum


# generate a vector of probabilities
def generate_alleles_probabilities(alleles_count):
    """
    Generate probabilities {p(i)}
    :param alleles_count: Amount of alleles
    :return: probabilities
    """
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
            if expected_val < cutoff:
                amount_of_small_expected_ += 1
                continue
            value += ((expected_val - observed_val) ** 2) / variance_val
    return value, amount_of_small_expected_


# generate simulated data as a list of lists (every row is: [k, i, j, p_k(i,j)])
def generate_data(alleles_count, population_amount, alleles_probabilities, alpha_val, uncertainty_val):
    """
    Generate a csv data simulation file
    :param alleles_count: Amount of alleles
    :param population_amount: Population size
    :param alleles_probabilities: Numpy array of alleles probabilities {p(i)}
    :param alpha_val: Float between 0.0 and 1.0, represents the closeness to HWE
    :param uncertainty_val: Float between 0.0 and 1.0, represents the uncertainty of the generated data
    :return: A csv file where every row represents: [id, first allele, second allele, probability]
    """

    # calculate probabilities {p(i,j)}
    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (
                        alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = 2 * alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
                # probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]

    # matrix 2 x N where every row contains the true alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, 2), dtype=np.int32)

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

    data = []
    # Vector of bernoulli's: for each person k, if choices[k] is 1 then this person will have uncertain alleles
    choices = random.choices(population=[0, 1], weights=[1 - uncertainty_val, uncertainty_val], k=population_amount)
    # Calculate the amount of donors with uncertain alleles
    sum_uncertain_donors = sum(choices)
    # The amount of pairs we will sample for each donor with uncertain alleles
    observations_amount_for_uncertain_donor = 10
    # Vector of size 10 * (amount of donors with uncertain alleles), every 10 values corresponds to an uncertain donor
    indices_uncertain = random.choices(population=range(len(probabilities_list)),
                                       weights=probabilities_list,
                                       k=observations_amount_for_uncertain_donor * sum_uncertain_donors)
    # This represents the starting index for the next 10 values in the vector defined above
    current_uncertainty_index = 0
    # Go over all the donors and add uncertainty to the uncertain donors
    for k in range(population_amount):
        # person k has true alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)
        # no uncertainty
        if choices[k] == 0:
            data.append([k, j, l, 1.0])
        else:
            # matrix 2 x 10 where each row represents 10 pairs of alleles for person k
            alleles_uncertainty = np.zeros(shape=(observations_amount_for_uncertain_donor, 2), dtype=int)
            # vector of probabilities for the uncertain paris of alleles for person k
            probabilities_uncertainty = np.zeros(shape=observations_amount_for_uncertain_donor)
            # Assign the uncertain alleles and probabilities
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
                # if t < 0.5 * observations_amount_for_uncertain_donor:
                #     probability *= 10000000
                # adding the weight of the uncertain alleles
                probabilities_uncertainty[t] = probability

            # normalizing the uncertain observations of current person into probabilities
            # probabilities_uncertainty = probabilities_uncertainty / np.sum(probabilities_uncertainty)

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
    difference_matrix = np.zeros(shape=(alleles_count, alleles_count))
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            mult = 1
            if i != j:
                mult = 2
            expected_val = mult * population_amount * alleles_probabilities[i] * alleles_probabilities[j]
            observed_val = population_amount * observed_probabilities[i, j]
            correction_val = correction[i, j]
            var_new = expected_val * correction_val
            var_old = expected_val

            i_j_new_stat = ((expected_val - observed_val) ** 2) / var_new
            i_j_old_stat = ((expected_val - observed_val) ** 2) / var_old
            if expected_val < cutoff_value:
                continue
            difference_matrix[i, j] = i_j_new_stat - i_j_old_stat
    difference_flat = difference_matrix.flatten()
    k_largest_indexes = (-difference_flat).argsort()[:20]
    for i in range(len(k_largest_indexes)):
        idx = k_largest_indexes[i]
        col = idx % difference_matrix.shape[1]
        row = (idx - col) // difference_matrix.shape[0]
        print(f'alleles: {row}, {col}, difference: {difference_matrix[row, col]}')

    dof_new = couples_amount - amount_of_small_expected_new

    p_value_new = 1 - stats.chi2.cdf(x=chi_squared_stat_new,
                                     df=dof_new)
    print(f'statistic old: {chi_squared_stat_old}, dof old: {dof_old}')
    print(f'statistic new: {chi_squared_stat_new}, dof new: {dof_new}')
    # return int(p_value_old < 0.05), int(p_value_new < 0.05), dof_old, dof_new
    return p_value_old, p_value_new, dof_old, dof_new


# given data: list of lists, run Traditional Chi Squared and Asta and return both significance result.
def plot_variance_vs_corrected_variance(data):
    id_to_index = {}
    allele_to_index = {}
    index_to_allele = {}
    # set of (i, j)
    ij_set = set()
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
        if (i, j) not in ij_set:
            ij_set.add((i, j))
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
    expected_corrected_vector = []
    variance_vector = []
    # for i in range(alleles_count):
    #     for j in range(i, alleles_count):
    #         expected = population_amount * alleles_probabilities[i] * alleles_probabilities[j]
    #         if i != j:
    #             expected *= 2
    #         expected_corrected = expected * correction[i, j]
    #         # if expected_corrected < 2.0:
    #         #     continue
    #         expected_corrected_vector.append(expected_corrected)
    index = 0
    for (i, j) in ij_set:
        print(f'loop: {index} / {len(ij_set)}')
        index += 1
        values_for_variance = []
        expected = population_amount * alleles_probabilities[i] * alleles_probabilities[j]
        if i != j:
            expected *= 2
        expected_corrected = expected * correction[i, j]
        # if expected_corrected < 2.0:
        #     continue
        expected_corrected_vector.append(expected_corrected)
        i_j = f'{i}_{j}'
        for current_id in id_to_i_j_to_i_j_observation:
            # for each id get the i_j observation, if not exists, get zero instead
            if i_j not in id_to_i_j_to_i_j_observation[current_id]:
                values_for_variance.append(0.0)
            else:
                values_for_variance.append(id_to_i_j_to_i_j_observation[current_id][i_j][2])
        mean_val = np.mean(values_for_variance)
        var_val = 0.0
        for element in values_for_variance:
            var_val += ((element - mean_val) ** 2)
        variance_vector.append(var_val)

    plt.scatter(expected_corrected_vector, variance_vector, color='deeppink')
    max_val = max(expected_corrected_vector + variance_vector)
    min_val = min(expected_corrected_vector + variance_vector)
    plt.xlim([min_val, max_val])
    plt.ylim([min_val, max_val])
    plt.xlabel('Expected * Corrected')
    plt.ylabel('Variance of (i,j) over {O_k(i,j)}_k')
    plt.show()

