import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D

if __name__ == '__main__':

    # fig, ax = plt.subplots(layout='constrained')
    fig, ax = plt.subplots()
    ax.grid(axis='y')
    ax.set_axisbelow(True)

    # arrays per test
    statistical_test_names = ['Traditional Chi-Squared',
                              'ASTA']
    # mean_paths = ['df_old_mean.csv',
    #               'df_correction_mean.csv',
    #               'df_sampling_mean.csv']
    mean_paths = ['means_old.csv',
                  'means_correction.csv']
    std_paths = ['stds_old.csv',
                 'stds_correction.csv']
    colors = ['blue',
              'red']

    # arrays per uncertainty value
    uncertainty_vals = [0.0, 0.1, 0.2, 0.4]
    markers = ['o',
               'D',
               's',
               'v']

    # alpha values
    interval_for_alpha = 0.04  # 0.02
    alpha_values = np.arange(start=0.92, stop=1 + interval_for_alpha,
                             step=interval_for_alpha)  # start, stop, step
    alphas = [str(alpha) for alpha in alpha_values]
    alphas = ['0.92', '0.96', '1.0']

    # list of lists (every list corresponds to a statistical test)
    transformations = []
    for row in range(len(statistical_test_names)):
        col = []
        for column in range(len(uncertainty_vals)):
            col.append(0)
        transformations.append(col)
    delta = 0.04
    offset = - delta * len(transformations) * len(transformations[0]) / 2
    for j in range(len(uncertainty_vals)):
        for i in range(len(statistical_test_names)):
            transformations[i][j] = Affine2D().translate(offset, 0.0) + ax.transData
            offset += delta

    for i in range(len(statistical_test_names)):
        for j, uncertainty in enumerate(uncertainty_vals):
            means_df = pd.read_csv(mean_paths[i], index_col=0)
            stds_df = pd.read_csv(std_paths[i], index_col=0)

            # get the jth column (according to uncertainty)
            means = means_df.iloc[:, j].tolist()
            stds = stds_df.iloc[:, j].tolist()

            # if uncertainty == 0.0:
            #     print(means_df)
            #     print(means)
            er1 = ax.errorbar(alphas, means, yerr=stds, marker=markers[j], color=colors[i], linestyle="none",
                              transform=transformations[i][j],
                              label=f'{statistical_test_names[i]}, uncertainty = {uncertainty}')
    plt.xlabel('Alpha values')
    plt.ylabel('Percentage of confidence amounts')
    plt.ylim(-0.2, 1.2)
    # plt.yticks(np.arange(0.0, 1.0, 0.1))
    plt.legend(bbox_to_anchor=(1.8, 1), loc='upper right', ncol=1)
    plt.savefig('testfig.png', bbox_inches='tight')
