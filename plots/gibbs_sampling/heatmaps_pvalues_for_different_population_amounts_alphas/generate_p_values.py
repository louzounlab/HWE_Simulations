import numpy as np
import pandas as pd
from models_simulated import gibbs_sampling

if __name__ == '__main__':
    population_amounts = np.array([100000, 1000000, 10000000, 100000000, 1000000000])
    # population_amounts = np.array([10000, 100000, 1000000, 10000000, 100000000])
    alpha_values = np.array([1.0, 0.99])

    alleles_amount = 100

    # initializing dataframes
    index = [f'{population_amount}' for population_amount in population_amounts]
    columns = [f'{alpha}' for alpha in alpha_values]

    df_p_values = pd.DataFrame(columns=columns, index=index, dtype=float)

    # filling dataframes
    for population_amount in population_amounts:
        for alpha in alpha_values:
            print(f'Current population amount: {population_amount}, current alpha: {alpha}')
            print('Preparing data')
            observed, _ \
                = gibbs_sampling.prepare_experiment_data(alleles_count=alleles_amount,
                                                         population_amount=population_amount,
                                                         alpha_val=alpha)
            print('Starting algorithm')
            p_value, elapsed_time = gibbs_sampling.full_algorithm(observations=observed)
            # we round the float to not have too many decimals
            df_p_values.loc[str(population_amount), str(alpha)] = p_value

    # saving dataframes
    df_p_values.to_csv('df_p_values.csv')

