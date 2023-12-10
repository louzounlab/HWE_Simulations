import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns
import utils


def plot_subfigure(index, row, col):
    ax = fig.add_subplot(gs[row, col])
    img = mpimg.imread(f'barplot_level_{levels_list[index]}.png')
    plt.title(levels_list[index])
    plt.imshow(img)


if __name__ == '__main__':
    sns.set()
    plt.rcParams["font.family"] = "Arial"

    levels_list = utils.LEVELS_LIST

    # fig = plt.figure(figsize=(10, 5))
    fig = plt.figure()

    rows_amount = 1
    cols_amount = 5
    gs = GridSpec(rows_amount, cols_amount, figure=fig)

    for i in range(len(levels_list)):
        col_ = i % cols_amount
        row_ = (i - col_) // rows_amount
        plot_subfigure(index=i, row=row_, col=col_)
    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=.0)
    plt.savefig('barplots_of_statistics_for_populations_races.svg', pad_inches=0.2, bbox_inches="tight", dpi=1000)
