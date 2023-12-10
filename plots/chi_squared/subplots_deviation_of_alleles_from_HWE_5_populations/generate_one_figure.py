import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns
import utils


def plot_subfigure(row, col):
    ax = fig.add_subplot(gs[row, col])
    img = mpimg.imread(f'plot_level_{levels_list[col]}_race_{races_list[row]}.png')
    plt.title(races_list[row])
    plt.imshow(img)


if __name__ == '__main__':
    sns.set()
    plt.rcParams["font.family"] = "Arial"

    levels_list = list(utils.LEVELS_LIST)
    races_list = list(utils.RACES_5_LIST)

    fig = plt.figure(figsize=(10, 10))

    gs = GridSpec(5, 5, figure=fig)

    for i in range(25):
        col_ = i % 5
        row_ = (i - col_) // 5
        plot_subfigure(row=row_, col=col_)
    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=.0)
    plt.savefig('deviation_of_alleles_from_HWE_5_populations.svg', pad_inches=0.2, bbox_inches="tight", dpi=1000)


