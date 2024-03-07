import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import seaborn as sns


if __name__ == '__main__':
    img = mpimg.imread('../heatmaps_pvalues_for_different_population_amounts_alphas/df_p_values.png')
    plt.imshow(img)
    plt.show()