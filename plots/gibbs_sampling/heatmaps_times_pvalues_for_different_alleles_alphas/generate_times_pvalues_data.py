import numpy as np
import pandas as pd
from models_simulated import gibbs_sampling

if __name__ == '__main__':
    alleles_nums = np.array([10, 20, 50, 100, 250, 500, 1000, 2000, 3000, 4000, 5000])
    alpha_values = np.array([1.0, 0.99])

    # LATEST RUN WITH ALPHA=0.99


    # alleles_nums = np.array([10])
    # alpha_values = np.array([1.0, 0.994])
    # 100M
    population_amount = 100000000
    # population_amount = 1000000

    # initializing dataframes
    index = [f'{allele}' for allele in alleles_nums]
    columns = [f'{alpha}' for alpha in alpha_values]

    df_times = pd.DataFrame(columns=columns, index=index, dtype=float)
    df_p_values = pd.DataFrame(columns=columns, index=index, dtype=float)

    # filling dataframes
    for alleles_amount in alleles_nums:
        for alpha in alpha_values:
            print(f'Current alleles amount: {alleles_amount}, current alpha: {alpha}')
            print('Preparing data')
            observed, _ \
                = gibbs_sampling.prepare_experiment_data(alleles_count=alleles_amount,
                                                         population_amount=population_amount,
                                                         alpha_val=alpha)
            print('Starting algorithm')
            p_value, elapsed_time = gibbs_sampling.full_algorithm(observations=observed)
            # we round the float to not have too many decimals
            df_times.loc[str(alleles_amount), str(alpha)] = elapsed_time
            df_p_values.loc[str(alleles_amount), str(alpha)] = p_value

    # saving dataframes
    df_times.to_csv('df_times.csv')
    df_p_values.to_csv('df_p_values.csv')