import matplotlib.pyplot as plt
import pandas as pd
import math
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns
import utils


def make_subplot(level_index, level):
    levels_list = utils.LEVELS_LIST
    races_5_list = utils.RACES_5_LIST
    races_21_list = utils.RACES_21_LIST
    races_1_list = utils.RACES_1_LIST

    index = range(21 + 5 + 1)
    columns = ['Population', 'Normalized Statistic']

    df_old = pd.DataFrame(index=index, columns=columns)
    df_new = pd.DataFrame(index=index, columns=columns)
    df_sampling = pd.DataFrame(index=index, columns=columns)
    counter = 0
    # Initialize the matplotlib figure
    # if level_index == 0:
    #     # in the first plot we have
    #     f, ax = plt.subplots(figsize=(6 * 502 / 602, 9))
    # else:
    #     f, ax = plt.subplots(figsize=(6, 9))
    for race in races_21_list:
        path_old_test = '../../../data/hla_data_21_populations/results/old_chi_squared_statistic'
        path_new_test = '../../../data/hla_data_21_populations/results/new_chi_squared_statistic'
        path_sampling_test = '../../../data/hla_data_21_populations/results/sampling_chi_squared_statistic'
        path_dof = '../../../data/hla_data_21_populations/results/dof'

        dof_data = pd.read_csv(path_dof, index_col=0)
        dof = dof_data.loc[level, race]

        statistics_data = pd.read_csv(path_old_test, index_col=0)
        old_statistic = statistics_data.loc[level, race] / dof

        statistics_data = pd.read_csv(path_new_test, index_col=0)
        new_statistic = statistics_data.loc[level, race] / dof

        statistics_data = pd.read_csv(path_sampling_test, index_col=0)
        sampling_statistic = statistics_data.loc[level, race] / dof

        df_old.iloc[counter] = [race, math.log(old_statistic, 10)]
        df_new.iloc[counter] = [race, math.log(new_statistic, 10)]
        df_sampling.iloc[counter] = [race, math.log(sampling_statistic, 10)]
        counter += 1

    for race in races_5_list:
        path_old_test = '../../../data/hla_data_5_populations/results/old_chi_squared_statistic'
        path_new_test = '../../../data/hla_data_5_populations/results/new_chi_squared_statistic'
        path_sampling_test = '../../../data/hla_data_5_populations/results/sampling_chi_squared_statistic'
        path_dof = '../../../data/hla_data_5_populations/results/dof'

        dof_data = pd.read_csv(path_dof, index_col=0)
        dof = dof_data.loc[level, race]

        statistics_data = pd.read_csv(path_old_test, index_col=0)
        old_statistic = statistics_data.loc[level, race] / dof

        statistics_data = pd.read_csv(path_new_test, index_col=0)
        new_statistic = statistics_data.loc[level, race] / dof

        statistics_data = pd.read_csv(path_sampling_test, index_col=0)
        sampling_statistic = statistics_data.loc[level, race] / dof

        df_old.iloc[counter] = [race, math.log(old_statistic, 10)]
        df_new.iloc[counter] = [race, math.log(new_statistic, 10)]
        df_sampling.iloc[counter] = [race, math.log(sampling_statistic, 10)]
        counter += 1

    for race in races_1_list:
        path_old_test = '../../../data/hla_data_1_populations/results/old_chi_squared_statistic'
        path_new_test = '../../../data/hla_data_1_populations/results/new_chi_squared_statistic'
        path_sampling_test = '../../../data/hla_data_1_populations/results/sampling_chi_squared_statistic'
        path_dof = '../../../data/hla_data_1_populations/results/dof'

        dof_data = pd.read_csv(path_dof, index_col=0)
        dof = dof_data.loc[level, race]

        statistics_data = pd.read_csv(path_old_test, index_col=0)
        old_statistic = statistics_data.loc[level, race] / dof

        statistics_data = pd.read_csv(path_new_test, index_col=0)
        new_statistic = statistics_data.loc[level, race] / dof

        statistics_data = pd.read_csv(path_sampling_test, index_col=0)
        sampling_statistic = statistics_data.loc[level, race] / dof

        df_old.iloc[counter] = [race, math.log(old_statistic, 10)]
        df_new.iloc[counter] = [race, math.log(new_statistic, 10)]
        df_sampling.iloc[counter] = [race, math.log(sampling_statistic, 10)]
        counter += 1

    # Plot the crashes where alcohol was involved
    # sns.set_color_codes("muted")
    # sns.barplot(x="Normalized Statistic", y="Population", data=df_new,
    #             label="new test", color="b")
    #
    # # Plot the total crashes
    # sns.set_color_codes("pastel")
    # sns.barplot(x="Normalized Statistic", y="Population", data=df_old,
    #             label="old test", color="b")
    df_new['ds'] = 'Asta'
    df_old['ds'] = 'Traditional Chi Squared'
    df_sampling['ds'] = 'Chi Squared with sampling'
    dss = pd.concat([df_old, df_new, df_sampling])
    # sns.barplot(x='Normalized Statistic', y='Population', hue='ds', data=dss)
    sns.barplot(x='Normalized Statistic', y='Population', hue='ds', data=dss)

    ax = plt.gca()
    # plt.xticks([])
    if level_index > 0:
        # ax.tick_params(left=False, bottom=False)
        plt.yticks([])
    else:
        # here we have the ticks for different races, we should color races from the same broad race
        races_to_broad_races = utils.races_to_broad_races
        broad_races_to_colors = utils.broad_races_to_colors

        colors = []
        for race in races_21_list:
            broad_race = races_to_broad_races[race]
            color = broad_races_to_colors[broad_race]
            colors.append(color)
        for broad_race in races_5_list:
            color = broad_races_to_colors[broad_race]
            colors.append(color)
        # we have the 'ALL' race left, only need to add a single color to it
        colors.append('gray')

        for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
            ticklabel.set_color(tickcolor)
    # ax.tick_params(axis='y', labelsize=10)

    # Add a legend and informative axis label
    plt.axvline(x=0.3, color='r', linestyle='--', label='Line: x = 0.3')
    ax.legend(fontsize=8)
    plt.savefig(f'barplot_level_{level}.pdf', format='pdf', bbox_inches="tight")


if __name__ == '__main__':
    sns.set_style('whitegrid')
    plt.rcParams["font.family"] = "Arial"

    levels_list = utils.LEVELS_LIST

    fig, axs = plt.subplots(1, 5, figsize=(16, 5), constrained_layout=True)

    for i, level_ in enumerate(levels_list):
        plt.subplot(1, 5, i + 1)
        plt.title(level_)
        make_subplot(level_index=i, level=level_)
    plt.savefig('barplots_of_statistics_for_populations_races.pdf', format='pdf', bbox_inches="tight")