import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import Affine2D
import seaborn as sns


def make_subplot(index, title, file_path):
    col = index % 2
    row = (index - col) // 2
    plt.subplot(gs[row, col])
    plt.title(title, weight='bold', fontsize=16)
    # plt.subplot(2, 2, index + 1)
    ax = plt.gca()
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
            means_df = pd.read_csv(f'{file_path}/{mean_paths[i]}', index_col=0)
            stds_df = pd.read_csv(f'{file_path}/{std_paths[i]}', index_col=0)

            # get the jth column (according to uncertainty)
            means = means_df.iloc[:, j].tolist()
            stds = stds_df.iloc[:, j].tolist()

            # if uncertainty == 0.0:
            #     print(means_df)
            #     print(means)
            er1 = ax.errorbar(alphas, means, yerr=stds, marker=markers[j], color=colors[i], linestyle="none",
                              transform=transformations[i][j],
                              label=f'{statistical_test_names[i]}, uncertainty = {uncertainty}')
    # plt.xlabel('Alpha values')
    # _, xlabels = plt.xticks()
    # _, ylabels = plt.yticks()
    # ax.set_xticklabels(xlabels, size=14)
    # ax.set_yticklabels(ylabels, size=14)
    plt.xlabel('Alpha values')
    plt.ylabel('Percentage of confidence amounts')
    # cax = ax.figure.axes[-1]
    # cax.tick_params(labelsize=14)
    # plt.ylabel('Percentage of confidence amounts')
    plt.ylim(-0.2, 1.2)
    # plt.yticks(np.arange(0.0, 1.0, 0.1))
    # plt.legend(bbox_to_anchor=(1.8, 1), loc='upper right', ncol=1)
    plt.legend()


if __name__ == '__main__':
    # sns.set()
    # fig, ax = plt.subplots()
    # fig.tight_layout()
    plt.rcParams["font.family"] = "Arial"
    # plt.figure(figsize=(9, 10))
    titles = ['A', 'B', 'C', 'D']
    alleles_amounts = [10, 15, 20, 25]  # 2
    population_sizes = [1000, 1000, 1000, 1000]  # 35

    file_paths = []
    for i in range(len(alleles_amounts)):
        alleles_amount = alleles_amounts[i]
        population_size = population_sizes[i]
        file_paths.append(f'population_{population_size}_alleles_{alleles_amount}')

    # first plot should take the first row
    # plt.subplot(2, 2, 1)
    fig = plt.figure(figsize=(14, 14), layout="constrained")
    gs = GridSpec(2, 2, figure=fig)

    # plot the other subplots
    for i in range(len(titles)):
        make_subplot(index=i, title=titles[i], file_path=file_paths[i])

    # plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    # adjust spacing
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=.15)
    plt.savefig('chi_squared_main_plot_supp_mat.svg', pad_inches=0.2, bbox_inches="tight", dpi=1000)
