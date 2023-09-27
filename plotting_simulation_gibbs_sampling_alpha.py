import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import single_experiment_with_certainty
import time

df_results_alpha_zero = pd.read_csv('data/simulation_data/gibbs_sampling_data/results_alpha=0', index_col=0)
df_results_alpha_fifty = pd.read_csv('data/simulation_data/gibbs_sampling_data/results_alpha=50', index_col=0)

plt.figure(figsize=(10, 8))
sns.heatmap(df_results_alpha_zero, annot=True, cmap=plt.cm.CMRmap_r)
plt.xlabel('Population')
plt.ylabel('Alleles')
plt.title('Means of % of higher results than 0. alpha = 0.0')
plt.show()

plt.figure(figsize=(10, 8))
sns.heatmap(df_results_alpha_fifty, annot=True, cmap=plt.cm.CMRmap_r)
plt.xlabel('Population')
plt.ylabel('Alleles')
plt.title('Means of % of higher results than 0. alpha = 0.5')
plt.show()