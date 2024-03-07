import numpy as np
import pandas as pd
from models_simulated import gibbs_sampling


# given alphas, alleles, population amounts: generate dataframes of p_values
def generate_p_values(alpha_values, alleles_nums, population_amounts):
    index = [f'{allele}' for allele in alleles_nums]
    columns = [f'{population_amount}' for population_amount in population_amounts]

    for i, alpha in enumerate(alpha_values):
        df_times = pd.DataFrame(columns=columns, index=index, dtype=float)

        for alleles_amount in alleles_nums:
            for population_amount in population_amounts:
                print(f'alleles: {alleles_amount}, population amount: {population_amount}')
                observed, _ \
                    = gibbs_sampling.prepare_experiment_data(alleles_count=alleles_amount,
                                                             population_amount=population_amount,
                                                             alpha_val=alpha)
                print('Starting algorithm')
                p_value, elapsed_time = gibbs_sampling.full_algorithm(observations=observed)
                # we round the float to not have too many decimals
                df_times.loc[str(alleles_amount), str(population_amount)] = elapsed_time

        df_times.to_csv(f'df_times.csv')


if __name__ == '__main__':
    # for large population
    alpha_values_ = np.array([1.0])
    alleles_nums_ = np.array([10, 20, 50, 100, 250, 500, 1000, 2000, 3000, 4000, 5000])
    population_amounts_ = np.array([1000000, 10000000, 100000000])
    generate_p_values(alpha_values=alpha_values_, alleles_nums=alleles_nums_,
                      population_amounts=population_amounts_)