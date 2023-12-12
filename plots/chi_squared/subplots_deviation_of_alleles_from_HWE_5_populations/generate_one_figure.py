import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns
import utils
from hwetests import asta


# def plot_subfigure(row, col):
#     ax = fig.add_subplot(gs[row, col])
#     img = mpimg.imread(f'plot_level_{levels_list[col]}_race_{races_list[row]}.png')
#     plt.title(races_list[row])
#     plt.imshow(img)


if __name__ == '__main__':
    sns.set_style('white')
    plt.rcParams["font.family"] = "Arial"

    levels_list = list(utils.LEVELS_LIST)
    races_list = list(utils.RACES_5_LIST)

    fig, axs = plt.subplots(5, 5, figsize=(8, 8), constrained_layout=True)

    counter = 1
    for level in levels_list:
        for race in races_list:
            plt.subplot(5, 5, counter)
            print(f'Current level: {level}, race: {race}')
            path = f'../../../data/hla_data_5_populations/levels/{level}/races/{race}/id_allele1_allele2_probability'
            file_to_save_name = f'plot_level_{level}_race_{race}.png'
            p_val, statistic, dof = asta.full_algorithm(file_path=path,
                                                        cutoff_value=2.0,
                                                        should_save_plot=file_to_save_name,
                                                        is_first_row_contains_columns_names=True,
                                                        title=race)
            counter += 1

    # for i in range(25):
    #     plt.subplot(5, 5, i + 1)
    #     col_ = i % 5
    #     row_ = (i - col_) // 5
    #     plot_subfigure(row=row_, col=col_)
    # plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=.0)
    plt.savefig('daviation_of_alleles_from_HWE_5_populations.pdf', format='pdf', bbox_inches="tight")


