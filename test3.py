import random
import numpy as np
import scipy.stats as stats
import os
import csv
import time
from models_simulated import chi_squared_second_attempt


def generate_data(alleles_count, population_amount, alleles_probabilities, alpha_val, uncertainty_val):
    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (
                        alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = 2 * alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
    print('probabilities matrix:')
    print(probabilities)
    print(f'sum of probabilities: {np.sum(probabilities)}')
    print('---------------------')

    alleles_individuals = np.zeros(
        (population_amount, 2), dtype=np.int32)

    # 1) assignment of alleles with certainty
    probabilities_list = probabilities.flatten()
    indices = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=population_amount)
    print(f'indices for certainty: {indices}')
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
    print('alleles_individuals:')
    print(alleles_individuals)
    print('---------------------')

    data = []
    choices = random.choices(population=[0, 1], weights=[1 - uncertainty_val, uncertainty_val], k=population_amount)
    print(f'uncertain choices: {choices}')
    sum_uncertain_donors = sum(choices)
    print(f'amount uncertain donors: {sum_uncertain_donors}')
    observations_amount_for_uncertain_donor = 1
    print(f'amount of pairs for each uncertain donor: {observations_amount_for_uncertain_donor}')
    indices_uncertain = random.choices(population=range(len(probabilities_list)),
                                       weights=probabilities_list,
                                       k=observations_amount_for_uncertain_donor * sum_uncertain_donors)
    print('indices for uncertainty:')
    print(indices_uncertain)
    print(f'amount of indices for uncertainty: {len(indices_uncertain)}')
    print('---------------------')

    current_uncertainty_index = 0

    # Go over all the donors and add uncertainty to the uncertain donors
    for k in range(population_amount):
        # person k has true alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)
        print(f'person {k} has certain alleles: {j}, {l}')
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
                probability = 1
                # if t < 0.5 * observations_amount_for_uncertain_donor:
                #     probability *= 10000000
                # adding the weight of the uncertain alleles
                probabilities_uncertainty[t] = probability

            probabilities_uncertainty = probabilities_uncertainty / np.sum(probabilities_uncertainty)
            for t in range(observations_amount_for_uncertain_donor):
                data.append([k,
                             alleles_uncertainty[t, 0],
                             alleles_uncertainty[t, 1],
                             probabilities_uncertainty[t]])
            current_uncertainty_index += observations_amount_for_uncertain_donor
    print('data:')
    print(data)
    return 1


alpha_value = 1.0
uncertainty_value = 0.4
population_size = 10
alleles_amount = 5

alleles_probability = chi_squared_second_attempt.generate_alleles_probabilities(alleles_count=alleles_amount)
print(f'alleles probability: {alleles_probability}')
print(f'sum is: {sum(alleles_probability)}')
print('---------------------')

data_ = generate_data(alleles_count=alleles_amount,
                      population_amount=population_size,
                      alleles_probabilities=alleles_probability,
                      alpha_val=alpha_value,
                      uncertainty_val=uncertainty_value)
