import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns


if __name__ == '__main__':
    sns.set()
    # fig, ax = plt.subplots()
    # fig.tight_layout()
    plt.rcParams["font.family"] = "Arial"

    titles = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    file_paths = ['hwe_design.png',
                  '../subplots_of_curves/subplots_of_curves_alpha_1000.png',
                  '../subplots_of_curves/subplots_of_curves_alpha_990.png',
                  '../heatmaps_times_pvalues_for_different_alleles_alphas/df_p_values.png',
                  '../heatmaps_pvalues_for_different_population_amounts_alphas/df_p_values.png',
                  '../scatter_times_for_different_population_amounts_alleles/df_times.png',
                  '../snipping_data/snipping_data_plot.png']

    fig = plt.figure(figsize=(9, 12), layout="constrained")

    gs = GridSpec(4, 2, figure=fig)

    ax0 = fig.add_subplot(gs[0, :])
    ax0.set_title(titles[0], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[0])
    plt.imshow(img)

    ax1 = fig.add_subplot(gs[1, 0])
    ax1.set_title(titles[1], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[1])
    ax1.imshow(img, aspect='auto')

    ax2 = fig.add_subplot(gs[1, 1])
    ax2.set_title(titles[2], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[2])
    ax2.imshow(img, aspect='auto')

    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title(titles[3], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[3])
    ax3.imshow(img, aspect='auto')

    ax4 = fig.add_subplot(gs[2, 1])
    ax4.set_title(titles[4], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[4])
    ax4.imshow(img, aspect='auto')

    ax5 = fig.add_subplot(gs[3, 0])
    ax5.set_title(titles[5], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[5])
    ax5.imshow(img, aspect='auto')

    ax6 = fig.add_subplot(gs[3, 1])
    ax6.set_title(titles[6], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[6])
    ax6.imshow(img, aspect='auto')

    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    # adjust spacing
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=.15)
    plt.savefig('gibbs_sampling_all_figures.svg', pad_inches=0.2, bbox_inches="tight", dpi=1000)
    plt.show()