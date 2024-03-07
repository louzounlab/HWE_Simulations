import random
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt


file_path = 'simulation_data.csv'
corrected_variances = []
observed_variances = []
with open(file_path, encoding="utf8") as infile:
    for index, line in enumerate(infile):
        # header
        if index == 0:
            continue
        lst = line.strip('\n').replace('+', ' ').replace(',', ' ').split()
        corrected_var = float(lst[5])
        observed_var = lst[6]
        if observed_var == '?':
            continue
        corrected_variances.append(corrected_var)
        observed_variances.append(float(observed_var))

list_values = corrected_variances + observed_variances
min_val = min(list_values)
max_value = max(list_values)

plt.scatter(corrected_variances, observed_variances, color='deeppink')
plt.xlabel('Corrected variance', fontsize=16)
plt.ylabel('Observed variance', fontsize=16)
plt.xlim([min_val, max_value])
plt.ylim([min_val, max_value])
ax = plt.gca()
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
plt.show()

