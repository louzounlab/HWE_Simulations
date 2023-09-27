import numpy as np
import random

import pandas as pd
import seaborn as sns
import array as arr

import statistical_test
import utils_with_certainty
import matplotlib.pyplot as plt
import time

alleles_count = 100
population_amount = 10
alpha_val = 0.0

alleles_probabilities = utils_with_certainty.calculate_alleles_probabilities(alleles_count)
print(alleles_probabilities)

observed = utils_with_certainty.calc_O_ij(population_amount, alleles_count, alpha_val, alleles_probabilities)

population_amount_calculated = 0
for i in range(alleles_count):
    for j in range(i, alleles_count):
        population_amount_calculated += observed[i, j]

# [ln(p_0), ln(p_1),...,]
# actually [1, delta_1, delta_2,...,delta_k]
list_probabilities = []

# calculate cdf. [(p_11, 1), (p_11+p_12, 2), ..., (p_11+...+p_kk, k)]
cdf_dict = utils_with_certainty.calculate_cdf_dict(alleles_count, observed)

# print(observed)
# print(cdf_dict)

# add first ln(probability) as 0. We only care about the deltas anyway.
list_probabilities.append(0)

# calc couples
first_couple = utils_with_certainty.calculate_couple_of_alleles(alleles_count, observed)
second_couple = random.choices(population=range(alleles_count), weights=alleles_probabilities,
                               k=2)