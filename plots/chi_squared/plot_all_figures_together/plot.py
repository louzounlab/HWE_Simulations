import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns


if __name__ == '__main__':
    sns.set()
    # fig, ax = plt.subplots()
    # fig.tight_layout()
    plt.rcParams["font.family"] = "Arial"

    titles = ['A', 'B', 'C', 'D', 'E', 'F']
    # file_paths = ['../main_plot/old_corrected_sampling.png',
    #               '../plot_var_no_sampling_vs_observed_and_var_corrected/variance_no_sampling_vs_observed_and_var_corrected_uncertainty_0.png',
    #               '../plot_var_no_sampling_vs_observed_and_var_corrected/variance_no_sampling_vs_observed_and_var_corrected_uncertainty_2.png',
    #               '../plot_var_no_sampling_vs_observed_and_var_corrected/variance_no_sampling_vs_observed_and_var_corrected_uncertainty_4.png',
    #               '../../testing_chi_squared_gibbs_with_sampling_ambiguity/df_p_values.png',
    #               '../deviation_of_one_allele_from_HWE_plot/deviation_of_one_allele_from_HWE_plot.png']

    file_paths = ['../../testing_chi_squared_gibbs_with_sampling_ambiguity/df_p_values.png',
                  '../plot_var_no_sampling_vs_observed_and_var_corrected/variance_no_sampling_vs_observed_and_var_corrected_uncertainty_0.png',
                  '../plot_var_no_sampling_vs_observed_and_var_corrected/variance_no_sampling_vs_observed_and_var_corrected_uncertainty_2.png',
                  '../plot_var_no_sampling_vs_observed_and_var_corrected/variance_no_sampling_vs_observed_and_var_corrected_uncertainty_4.png',
                  '../main_plot/main_plot.png',
                  '../deviation_of_one_allele_from_HWE_plot/deviation_of_one_allele_from_HWE_plot.png']

    # fig = plt.figure(figsize=(7, 8.8))
    fig = plt.figure(figsize=(5, 2.5), layout="constrained")
    # plt.subplot(gs[0, :2])
    gs = GridSpec(2, 4, figure=fig)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.set_title(titles[0], weight='bold', fontsize=12)
    img = mpimg.imread(file_paths[0])
    plt.imshow(img)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.set_title(titles[1], weight='bold', fontsize=12)
    img = mpimg.imread(file_paths[1])
    ax1.imshow(img)

    ax2 = fig.add_subplot(gs[0, 2])
    ax2.set_title(titles[2], weight='bold', fontsize=10)
    img = mpimg.imread(file_paths[2])
    ax2.imshow(img)

    ax3 = fig.add_subplot(gs[0, 3])
    ax3.set_title(titles[3], weight='bold', fontsize=10)
    img = mpimg.imread(file_paths[3])
    ax3.imshow(img)

    ax4 = fig.add_subplot(gs[1, :3])
    ax4.set_title(titles[4], weight='bold', fontsize=10)
    img = mpimg.imread(file_paths[4])
    ax4.imshow(img)

    ax5 = fig.add_subplot(gs[1, 3])
    ax5.set_title(titles[5], weight='bold', fontsize=10)
    img = mpimg.imread(file_paths[5])
    ax5.imshow(img)

    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    # adjust spacing
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=.0)
    plt.savefig('chi_squared_all_figures.svg', pad_inches=0.2, bbox_inches="tight", dpi=1000)
    # plt.show()