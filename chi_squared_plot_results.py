import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


names = ['ambiguities', 'new_chi_squared_p_values', 'new_chi_squared_statistic', 'old_chi_squared_p_values',
         'old_chi_squared_statistic', 'dof']

for name in names:
    df = pd.read_csv(f'data/real_data_test/results/{name}', index_col=0)

    plt.figure(figsize=(10, 8))
    sns.heatmap(df, annot=True, cmap=plt.cm.CMRmap_r)
    plt.xlabel('Race')
    plt.ylabel('Locus')
    plt.title(name)
    plt.show()