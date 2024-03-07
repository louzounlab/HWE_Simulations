import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns


if __name__ == '__main__':
    sns.set()
    # fig, ax = plt.subplots()
    # fig.tight_layout()
    plt.rcParams["font.family"] = "Arial"

    titles = ['A', 'B']
    file_paths = ['df_p_values_small_population_alpha_920.png',
                  'df_p_values_small_population_alpha_1000.png']

    fig = plt.figure(layout="constrained")

    gs = GridSpec(1, 2, figure=fig)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.set_title(titles[0], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[0])
    plt.imshow(img)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.set_title(titles[1], weight='bold', fontsize=14)
    img = mpimg.imread(file_paths[1])
    ax1.imshow(img)

    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])

    plt.savefig('gibbs_sampling_small_population.svg', pad_inches=0.2, bbox_inches="tight", dpi=1000)
    plt.show()