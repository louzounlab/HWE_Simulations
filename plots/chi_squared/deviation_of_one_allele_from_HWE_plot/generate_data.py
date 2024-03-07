import numpy as np
from matplotlib import pyplot as plt
import random
import csv


# assuming x is 2-dimensional
def softmax_1d(x):
    sum_ = 0.0
    for row in range(x.shape[0]):
        sum_ += np.exp(x[row])
    return np.exp(x) / sum_


# generate a vector of probabilities
def calculate_alleles_probabilities(alleles_count):
    probs = np.random.uniform(0.0, 2.0, size=alleles_count)
    probs = softmax_1d(probs)
    return probs


def prepare_probabilities(alleles_count):
    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)
    return alleles_probabilities


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   alleles_probabilities):
    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    # the allele we choose to be outside of hwe
    allele_outside_hwe = 0
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = alleles_probabilities[t] * alleles_probabilities[m]
                probabilities[m, t] = alleles_probabilities[t] * alleles_probabilities[m]

    # add deviation from HWE for the chosen allele (here m ranges over all the alleles)
    for m in range(alleles_count):
        y, u = min(m, allele_outside_hwe), max(m, allele_outside_hwe)
        if y == u:
            probabilities[y, y] = (2 * alpha_val - 1) * (alleles_probabilities[y] ** 2) + 2 * (1 - alpha_val) * alleles_probabilities[y]
        else:
            # we don't multiply by 2 here yet
            probabilities[y, u] = alpha_val * alleles_probabilities[y] * alleles_probabilities[u]
            probabilities[u, y] = alpha_val * alleles_probabilities[y] * alleles_probabilities[u]
    print(f'sum of probs: {np.sum(probabilities)}')
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

    # create csv
    file_name = 'data.csv'
    columns = ['id', 'first_allele', 'second_allele', 'probability']

    with open(file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(columns)

        # adding uncertainty to our model
        choices = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=population_amount)
        for k in range(population_amount):
            if k % 1000 == 0:
                print(f'person: {k} / {population_amount}')
            # person k has alleles j,l
            j, l = alleles_individuals[k]
            j, l = min(j, l), max(j, l)

            # choice whether this person will have uncertain alleles
            choice = choices[k]

            # this person has uncertain alleles
            if choice == 1:
                writer.writerow([k, j, l, 1.0])
            # this person has certain alleles
            if choice == 0:
                for y in range(alleles_count):
                    for u in range(y, alleles_count):
                        writer.writerow([k, y, u, probabilities[y, u]])


if __name__ == '__main__':
    alleles_amount = 100
    population = 10000
    alpha = 0.7
    uncertainty = 0.2
    alleles_probabilities = calculate_alleles_probabilities(alleles_amount)

    run_experiment(alleles_count=alleles_amount,
                   population_amount=population,
                   alpha_val=alpha,
                   uncertainty_val=uncertainty,
                   alleles_probabilities=alleles_probabilities)

