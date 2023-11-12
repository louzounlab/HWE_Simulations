import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns


def make_subplot(index, title, file_path):
    plt.subplot(4, 2, index + 1)
    plt.title(title, weight='bold', fontsize=14)
    img = mpimg.imread(file_path)
    plt.imshow(img, aspect='auto')


if __name__ == '__main__':
    sns.set()
    # fig, ax = plt.subplots()
    # fig.tight_layout()
    plt.rcParams["font.family"] = "Arial"
    plt.figure(figsize=(9, 10))

    titles = ['A', 'B', 'C', 'D', 'E', 'F']
    file_paths = ['hwe_design.png',
                  '../subplots_of_curves/subplots_of_curves_alpha_1000.png',
                  '../subplots_of_curves/subplots_of_curves_alpha_990.png',
                  '../heatmaps_times_pvalues_for_different_alleles_alphas/df_p_values.png',
                  '../heatmaps_pvalues_for_different_population_amounts_alleles/df_p_values_small_population_alpha_920.png',
                  '../scatter_times_for_different_population_amounts_alleles/df_times.png',
                  '../snipping_data/snipping_data_plot.png']

    # first plot should take the first row
    plt.subplot(1, 2, 1)
    plt.title(titles[0], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[0])
    plt.imshow(img, aspect='auto')

    # plot the other subplots
    for i in range(1, len(titles)):
        make_subplot(index=i, title=titles[i], file_path=file_paths[i])

    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    # adjust spacing
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=.15)
    plt.savefig('image.svg', pad_inches=0.2, bbox_inches="tight", dpi=1000)
    plt.show()
