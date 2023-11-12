import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
from models_simulated import gibbs_sampling


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


def prepare_probabilities(alleles_count):
    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)

    return alleles_probabilities


def get_O_ij(observed, i, j):
    row, col = min(i, j), max(i, j)
    return observed[row, col]


def get_p_ij(probabilities, i, j):
    i, j = min(i, j), max(i, j)
    return probabilities[i, j]


def calc_observed_cdf(alleles_count, observed):
    # list of numpy arrays
    observed_cdf = []
    for i in range(alleles_count):
        np_array = np.zeros(alleles_count)
        # 1,i 2,i 3,i ,..., i,i i,i+1 i,i+2 ,..., i,n
        for k in range(alleles_count):
            t, m = min(i, k), max(i, k)
            mult = 0.5
            if t == m:
                mult = 1.0
            np_array[k] = mult * observed[t, m]
        observed_cdf.append(np_array)
    return observed_cdf


# calculate the starting index for the gibbs sampling
def calculate_start_time(alleles_count, population_amount, alleles_probabilities, observed, observed_cdf, iterations):
    observed_cdf = [np_array / sum(np_array) for np_array in observed_cdf]
    start_time = 0.0
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            mult = 1
            if i != j:
                mult = 2
            val = observed[i, j] - mult * population_amount * alleles_probabilities[i] * alleles_probabilities[j]
            if val >= 0:
                val = abs(val) * observed_cdf[i][j] * (1 - alleles_probabilities[j])
            else:
                val = abs(val) * (1 - observed_cdf[i][j]) * alleles_probabilities[j]
            start_time += abs(val)
    start_time = int(start_time)
    if start_time > 50000:
        return 50000
    if start_time < 10000:
        return 10000
    return start_time


def update_current_delta_probability(list_probabilities: list, observed, observed_cdf, alleles_probabilities, couples,
                                     population_amount_calculated, probabilities):
    for row in range(couples.shape[0]):
        # here we make sure that i <= j so we can access the upper triangle of probabilities
        i, j = min(couples[row, 0], couples[row, 1]), max(couples[row, 0], couples[row, 1])
        # i, j = int(i), int(j)

        # z = np.log(population_amount_calculated * probabilities[i, j] / (1 - probabilities[i, j]))
        sign = couples[row, 2]
        val = 0
        try:
            if sign == -1:
                # val = np.log(observed[i, j] / (population_amount_calculated - observed[i, j] + 1)) \
                #       + np.log((1 - get_p_ij(alleles_probabilities, i, j)) / get_p_ij(alleles_probabilities, i, j))
                val = np.log(get_O_ij(observed, i, j)) - np.log(
                    population_amount_calculated - get_O_ij(observed, i, j) + 1) \
                      - np.log(get_p_ij(probabilities, i, j)) + np.log(1 - get_p_ij(probabilities, i, j))
            elif sign == 1:
                val = np.log(population_amount_calculated - get_O_ij(observed, i, j)) - np.log(
                    get_O_ij(observed, i, j) + 1) \
                      + np.log(get_p_ij(probabilities, i, j)) - np.log(1 - get_p_ij(probabilities, i, j))
        except RuntimeWarning:
            print(f'''
i,j: {i}, {j}.
O_ij: {get_O_ij(observed, i, j)}
population calculated: {population_amount_calculated}
p_ij: {get_p_ij(probabilities, i, j)}
sign: {sign}
observations: {observed}
alleles probabilities: {alleles_probabilities}''')
            raise ValueError('A very specific bad thing happened.')

        # print(observed)
        # print(f'i: {i}, j: {j}')
        # sum_current += sign * (z + np.log(observed[i, j]))
        list_probabilities.append(val)
        # unbalanced_alleles[i] += sign
        # unbalanced_alleles[j] += sign
        observed[i, j] += sign

        # also update observed_cdf
        # o(i,j)+1 -> j_th row, i_th row
        if i == j:
            observed_cdf[j][i] += sign
        else:
            observed_cdf[j][i] += 0.5 * sign
            observed_cdf[i][j] += 0.5 * sign


# given the simulated data (or real data), perform Gibbs Sampling.
# returns:  result (1 for significance or 0)
def perform_experiment(alleles_count,
                       population_amount_calculated, alleles_probabilities, probabilities, observed, observed_cdf,
                       plot_index=-1):
    # [ln(p_0), ln(p_1),...,]
    # actually [0, delta_1, delta_2,...,delta_k]
    list_probabilities = []

    range_alleles_count = range(alleles_count)

    # print(observed)

    # calculate cdf. [(p_11, 1), (p_11+p_12, 2), ..., (p_11+...+p_kk, k)]
    # cdf_dict = {}
    # utils_with_certainty.calculate_cdf_dict(alleles_count, observed, cdf_dict=cdf_dict)

    # print(cdf_dict)

    # add first ln(probability) as 0. We only care about the deltas anyway.
    list_probabilities.append(0)

    # calc couples
    # first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
    # first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
    i = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                       k=1)[0]

    j = random.choices(population=range_alleles_count, weights=observed_cdf[i],
                       k=1)[0]
    first_couple = [i, j]

    k = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                       k=1)[0]
    # while k in {i, j}:
    #     k = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                        k=1)[0]
    l = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                       k=1)[0]
    # l = k
    # while l in {i, j}:
    #     l = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                        k=1)[0]
    second_couple = [k, l]

    # every row represents a couple i,j and a value (either 1 or -1)
    # in the first iteration here we have 4 couples, but in the next iterations we will have 2 couples.
    couples = np.zeros(shape=(2, 3), dtype=int)
    couples[0, :] = [first_couple[0], first_couple[1], -1]
    couples[1, :] = [second_couple[0], second_couple[1], +1]

    # calculate new delta and modify observed
    update_current_delta_probability(list_probabilities, observed, observed_cdf, alleles_probabilities, couples,
                                     population_amount_calculated,
                                     probabilities)

    # pick two alleles using the cdf (first we pick a number between 0 and 1 and then get element with the closest
    # probability)
    # t
    allele_1 = random.choices(population=range_alleles_count, weights=alleles_probabilities,
                              k=1)[0]
    # while allele_1 in {first_couple[0], first_couple[1], second_couple[1]}:
    #     allele_1 = random.choices(population=range(alleles_count), weights=alleles_probabilities,
    #                               k=1)[0]
    # m
    allele_2 = random.choices(population=range_alleles_count, weights=observed_cdf[second_couple[1]],
                              k=1)[0]
    # while allele_2 in {allele_1}:
    #     allele_2 = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
    #                                                              left_alleles=second_couple)

    couple_from_cdf = [allele_1, allele_2]

    # every row represents a couple i,j and a value (either 1 or -1)
    # in the first iteration here we have 4 couples, but in the next iterations we will have 2 couples.
    couples[0, :] = [first_couple[1], allele_1, +1]
    couples[1, :] = [second_couple[1], allele_2, -1]

    # calculate new delta and modify observed
    update_current_delta_probability(list_probabilities, observed, observed_cdf, alleles_probabilities, couples,
                                     population_amount_calculated,
                                     probabilities)
    # now keep iterating to fill the list of deltas
    # iterations: 100000
    iterations = 100000
    # iterations = 5
    for k in range(iterations):
        # sums.append(sum_loops[0])
        # if (k % 10000) == 0:
        #     print(f' loop {k} / {iterations}. experiment num: . alleles: {alleles_count}. '
        #           f'population: {population_amount}. alpha: ')
        # print(observed)
        # print('cdf dict:')
        # print(cdf_dict)
        # print(alleles_state)

        # updating the couples
        first_couple = [first_couple[0], couple_from_cdf[1]]
        second_couple = [second_couple[0], couple_from_cdf[0]]

        # pick two alleles using the cdf (first we pick a number between 0 and 1 and then get element with the closest
        # probability)
        # x
        x = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                           k=1)[0]
        # while x in {first_couple[0], first_couple[1], second_couple[1]}:
        #     x = random.choices(population=range(alleles_count), weights=alleles_probabilities,
        #                        k=1)[0]
        # y
        y = random.choices(population=range_alleles_count, weights=observed_cdf[second_couple[1]],
                           k=1)[0]
        # while y in {x}:
        #     y = utils_with_certainty.get_allele_from_cdf_dict(alleles_count, cdf_dict, observed,
        #                                                       left_alleles=second_couple)

        couple_from_cdf = [x, y]

        couples[0, :] = [first_couple[1], x, +1]
        couples[1, :] = [second_couple[1], y, -1]
        # print(couples)

        # print(f'plus couple to update: {first_couple[1]}, {x}')
        # print(f'minus couple to update: {second_couple[1]}, {y}')

        update_current_delta_probability(list_probabilities, observed, observed_cdf, alleles_probabilities,
                                         couples,
                                         population_amount_calculated,
                                         probabilities)
        # print(observed)

    # now we have the list of probabilities. check if 95% of the elements (sum of deltas) are bigger than 1.
    sum_current = 0
    bigger_counter = 0
    # start_from = 10000
    start_from = calculate_start_time(alleles_count=alleles_count,
                                      population_amount=population_amount_calculated,
                                      alleles_probabilities=alleles_probabilities,
                                      observed=observed,
                                      observed_cdf=observed_cdf,
                                      iterations=iterations)

    values = []
    for delta in list_probabilities[:start_from]:
        sum_current += delta

    for delta in list_probabilities[start_from:]:
        sum_current += delta
        values.append(sum_current)
        if sum_current >= list_probabilities[0]:
            bigger_counter += 1

        # if plot_index is zero we have only one plot, if bigger than zero we have subplots.
        if plot_index > 0:
            plt.subplot(3, 2, plot_index)
            plt.xticks([])
            if (plot_index % 2) == 0:
                plt.yticks([])
            ax = plt.gca()
            ax.tick_params(axis='y', labelsize=10)

        if plot_index >= 0:
            plt.plot(list(range(start_from, len(values) + start_from)), values, color='black')

    result = bigger_counter / (len(list_probabilities) - start_from)
    p_value = 1.0 - 2 * abs(max(0.5, result) - 0.5)
    return p_value


# here is the full gibbs sampling algorithm.
# we use the same observations 5 times and return the mean result
def full_algorithm(observations, should_save_plot=False):
    alleles_count = observations.shape[0]
    population_amount_calculated = np.sum(observations)
    alleles_probabilities = np.zeros(alleles_count)

    for i in range(alleles_count):
        probability = 0.0
        for j in range(alleles_count):
            if i == j:
                probability += observations[i, i]
            else:
                probability += observations[min(i, j), max(i, j)] * 0.5
        alleles_probabilities[i] = probability / population_amount_calculated

    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            mult = 2.0
            if i == j:
                mult = 1.0
            probabilities[i, j] = mult * alleles_probabilities[i] * alleles_probabilities[j]

    observed_cdf = calc_observed_cdf(alleles_count, observations)

    experiments_amount = 1
    observed_copy = np.copy(observations)
    observed_cdf_copy = np.copy(observed_cdf)
    results = []

    for experiment_num in range(experiments_amount):
        if should_save_plot:
            plt.subplot(3, 2, 1)
            result = \
                perform_experiment(alleles_count=alleles_count,
                                   population_amount_calculated=population_amount_calculated,
                                   alleles_probabilities=alleles_probabilities,
                                   probabilities=probabilities,
                                   observed=observed_copy,
                                   observed_cdf=observed_cdf_copy,
                                   plot_index=experiment_num + 1)
        else:
            result = \
                perform_experiment(alleles_count=alleles_count,
                                   population_amount_calculated=population_amount_calculated,
                                   alleles_probabilities=alleles_probabilities,
                                   probabilities=probabilities,
                                   observed=observed_copy,
                                   observed_cdf=observed_cdf_copy)
        results.append(result)
        # initialize observed_copy, copy from observed matrix
        np.copyto(observed_copy, observations)
        np.copyto(observed_cdf_copy, observed_cdf)

    if should_save_plot:
        if isinstance(should_save_plot, str):
            plt.savefig(should_save_plot, pad_inches=0.2, bbox_inches="tight")
        else:
            plt.savefig(f'gibbs_sampling_plot', pad_inches=0.2, bbox_inches="tight")
    mean_result = np.mean(results)
    return mean_result


def prepare_data(alleles_count, population_amount, alpha_val, uncertainty_val,
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
    all_probabilities = np.zeros(shape=(alleles_count, alleles_count, population_amount))

    choices = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=population_amount)

    # adding uncertainty to our model
    for k in range(population_amount):
        # person k has alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)

        # choice whether this person will have uncertain alleles
        choice = choices[k]

        # this person has certain alleles
        if choice == 1:
            all_probabilities[j, l, k] = 1.0
        # this person has uncertain alleles
        if choice == 0:
            all_probabilities[:, :, k] = probabilities

    # for each k make {p_k(i,j)} symmetric and choose for every k alleles i,j
    for k in range(population_amount):
        for t in range(alleles_count):
            for m in range(t, alleles_count):
                all_probabilities[t, m, k] = (all_probabilities[t, m, k] + all_probabilities[m, t, k]) / 2
                all_probabilities[m, t, k] = (all_probabilities[t, m, k] + all_probabilities[m, t, k]) / 2
        probabilities_list = all_probabilities[:, :, k].flatten()
        index = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=1)[0]

        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]
    # calculate observations
    observations = np.zeros(shape=(alleles_count, alleles_count))
    for k in range(population_amount):
        t, m = alleles_individuals[k]
        observations[t, m] += 1
    return observations


if __name__ == '__main__':
    alleles_amount = 100
    alpha_values = np.array([1.0, 0.92])
    uncertainty_values = np.array([0.0, 0.2, 0.4])

    # LATEST RUN WITH ALPHA=0.99

    # alleles_nums = np.array([10])
    # alpha_values = np.array([1.0, 0.994])
    # 100M
    population_amount = 10000
    # population_amount = 1000000

    # initializing dataframes
    index = [f'{uncertainty}' for uncertainty in uncertainty_values]
    columns = [f'{alpha}' for alpha in alpha_values]

    df_p_values = pd.DataFrame(columns=columns, index=index, dtype=float)

    alleles_probabilities_ = prepare_probabilities(alleles_count=alleles_amount)

    # filling dataframes
    for alpha in alpha_values:
        for uncertainty in uncertainty_values:
            print(f'Current alpha: {alpha}, current uncertainty: {uncertainty}')

            observations_ = prepare_data(alleles_count=alleles_amount, population_amount=population_amount,
                                         alpha_val=alpha,
                                         uncertainty_val=uncertainty,
                                         alleles_probabilities=alleles_probabilities_)
            print('Starting algorithm')
            p_value, elapsed_time = gibbs_sampling.full_algorithm(observations=observations_)
            df_p_values.loc[str(uncertainty), str(alpha)] = p_value

    # saving dataframes
    df_p_values.to_csv('df_p_values.csv')
