import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import single_experiment_with_certainty
import time

if __name__ == '__main__':
    np.seterr(all='raise')
    alleles_nums = np.array([10])
    population_nums = np.array([100])
    experiments_amount = 1
    alpha_vals = [0.0]

    for alpha_val in alpha_vals:
        for i, alleles in enumerate(alleles_nums):
            for j, population in enumerate(population_nums):
                # prepare the simulation data
                population_amount_calculated, alleles_probabilities, probabilities, observed, elapsed_initial_time = \
                    single_experiment_with_certainty.prepare_experiment_data(alleles_count=alleles,
                                                                             population_amount=population,
                                                                             alpha_val=alpha_val)

                single_experiment_with_certainty.perform_experiment(alleles_count=alleles,
                                                                    population_amount=population,
                                                                    alpha_val=alpha_val,
                                                                    experiment_num=experiments_amount,
                                                                    population_amount_calculated=population_amount_calculated,
                                                                    alleles_probabilities=alleles_probabilities,
                                                                    probabilities=probabilities,
                                                                    observed=observed)