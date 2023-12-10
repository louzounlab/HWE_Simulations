import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from models_simulated import chi_squared

fix, ax = plt.subplots(figsize=(17, 6))
ax.grid()
ax.set_axisbelow(True)

warnings.filterwarnings("error")

num_of_experiments = 1000  # 500 amount of experiments for each alpha ranging from zero to one
alleles_amount = 10  # 2
population_size = 2000  # 35
interval_for_alpha = 0.04  # 0.02
# uncertainty = 0.2
uncertainty_vals = [0.0, 0.1, 0.2, 0.4]
# blue color: 'royalblue'
colors = ['royalblue', 'slategrey', 'limegreen', 'deeppink']
alpha_values = np.arange(start=0.92, stop=1 + interval_for_alpha,
                         step=interval_for_alpha)  # start, stop, step
alpha_values = np.array([round(alpha, 2) for alpha in alpha_values])
# alpha_values = np.array([1.0])
confidence_amounts_per_alpha_old_test = np.zeros((len(uncertainty_vals), alpha_values.shape[0]))
confidence_amounts_per_alpha_correction_test = np.zeros((len(uncertainty_vals), alpha_values.shape[0]))
confidence_amounts_per_alpha_sampling_test = np.zeros((len(uncertainty_vals), alpha_values.shape[0]))

df_old = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
df_correction = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
df_sampling = pd.DataFrame(index=alpha_values, columns=uncertainty_vals)
should_use_new_test = 0
# markers = ['', '>']
markers = ['.', '>', 'o']
labels = ['Old Chi-Squared', 'Chi Squared Correction', 'Chi Squared Sampling']
# print(alpha_values.shape)
# print(confidence_amounts_per_alpha.shape)
alleles_probabilities_ = chi_squared.prepare_probabilities(alleles_count=alleles_amount)

for uncertainty_index, uncertainty in enumerate(uncertainty_vals):
    for i, alpha in enumerate(alpha_values):
        counter_old_test = 0
        counter_correction_test = 0
        counter_sampling_test = 0

        for experiment in range(num_of_experiments):
            result_old, result_correction, result_sampling = chi_squared.run_experiment(alleles_amount, population_size,
                                                                                        alpha, uncertainty,
                                                                                        alleles_probabilities_)
            counter_old_test += result_old
            counter_correction_test += result_correction
            counter_sampling_test += result_sampling
        print(f'uncertainty = {uncertainty}, alpha = {alpha}')
        confidence_amount_old = counter_old_test / num_of_experiments
        confidence_amount_correction = counter_correction_test / num_of_experiments
        confidence_amount_sampling = counter_sampling_test / num_of_experiments

        confidence_amounts_per_alpha_old_test[uncertainty_index, i] = confidence_amount_old
        df_old.iloc[i, uncertainty_index] = confidence_amount_old

        confidence_amounts_per_alpha_correction_test[uncertainty_index, i] = confidence_amount_correction
        df_correction.iloc[i, uncertainty_index] = confidence_amount_correction

        confidence_amounts_per_alpha_sampling_test[uncertainty_index, i] = confidence_amount_sampling
        df_sampling.iloc[i, uncertainty_index] = confidence_amount_sampling

    # print(confidence_amounts_per_alpha[uncertainty_index, :])
    plt.plot(alpha_values, confidence_amounts_per_alpha_old_test[uncertainty_index, :],
             label=f'{labels[0]}, uncertainty = {uncertainty}',
             marker=markers[0],
             color=colors[uncertainty_index])

    plt.plot(alpha_values, confidence_amounts_per_alpha_correction_test[uncertainty_index, :],
             label=f'{labels[1]}, uncertainty = {uncertainty}',
             marker=markers[1],
             color=colors[uncertainty_index])

    plt.plot(alpha_values, confidence_amounts_per_alpha_sampling_test[uncertainty_index, :],
             label=f'{labels[2]}, uncertainty = {uncertainty}',
             marker=markers[2],
             color=colors[uncertainty_index])

# print(confidence_amounts_per_alpha[-1])
df_old.to_csv('df_old.csv')
df_correction.to_csv('df_correction.csv')
df_sampling.to_csv('df_sampling.csv')

# plt.plot(alpha_values, confidence_amounts_per_alpha, 'r-', label="")
plt.xlabel('Alpha values', fontsize=14)
plt.ylabel('Percentage of confidence amounts', fontsize=14)
plt.yticks(np.arange(0.0, 1.0, 0.05))
ax = plt.gca()
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
# plt.title('Old vs Corrected Chi-Squared', fontsize=14)
plt.legend()
# plot, plot_opposite, plot_small
plt.savefig('old_corrected_sampling', pad_inches=0.2, bbox_inches="tight")
plt.show()
