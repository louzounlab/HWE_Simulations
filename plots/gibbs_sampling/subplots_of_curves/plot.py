import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def make_plot(alpha_str):

    for i in range(1, 50):
        df = pd.read_csv(f'alpha_{alpha_str}_plotindex_{i}', index_col=0)
        times = df.iloc[:, 0].to_list()
        values = df.iloc[:, 1].to_list()
        plt.subplot(7, 7, i)
        plt.xticks([])
        if i % 7 != 1:
            plt.yticks([])
        ax = plt.gca()
        ax.tick_params(axis='y', labelsize=12)
        plt.plot(times, values, color='black')


if __name__ == '__main__':

    alpha_vals = [1.0, 0.99]

    for alpha_val in alpha_vals:
        sns.set_style('white')
        # sns.set(font_scale=1.3)
        plt.rcParams["font.family"] = "Arial"
        fig, axs = plt.subplots(figsize=(9, 6))
        alpha_str_val = str(int(alpha_val * 100))
        make_plot(alpha_str_val)
        plt.savefig(f'subplots_of_curves_alpha_{alpha_str_val}.pdf', format='pdf', bbox_inches="tight")