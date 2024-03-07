import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
from models_simulated import gibbs_sampling
from models_simulated import chi_squared_second_attempt


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
