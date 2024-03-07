import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def softmax_1d(x):
    sum = 0.0
    for row in range(x.shape[0]):
        sum += np.exp(x[row])
    return np.exp(x) / sum


# generate a vector of probabilities
def calculate_alleles_probabilities(alleles_count):
    probs = np.random.uniform(0.0, 2.0, size=alleles_count)
    probs = softmax_1d(probs)
    return probs


def calc_variance(input_list, input_mean):
    var = 0.0
    for element in input_list:
        var += (element - input_mean) ** 2
    return var


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   alleles_probabilities, plot_index):
    probabilities = np.zeros(shape=(alleles_count, alleles_count))
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            if t == m:
                probabilities[t, m] = (1 - alpha_val) * alleles_probabilities[t] + alpha_val * (
                        alleles_probabilities[m] ** 2)
            else:
                # we don't multiply by 2 here yet
                probabilities[t, m] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]
                probabilities[m, t] = alpha_val * alleles_probabilities[t] * alleles_probabilities[m]

    # matrix where every row is the alleles for a person
    alleles_individuals = np.zeros(
        (population_amount, 2), dtype=np.int32)

    # print(f'alleles: {alleles_probabilities}')
    # print(f' marginal: {marginal_probabilities}')

    # 1) assignment of alleles with certainty
    probabilities_list = probabilities.flatten()
    indices = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=population_amount)
    for k in range(population_amount):
        # notice the alleles i != j have probability 2 * p(i,j)
        index = indices[k]
        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

        # matrix p_k(i,j) = A[i,j,k]
    all_probabilities = np.zeros(shape=(alleles_count, alleles_count, population_amount))

    choices = random.choices(population=[0, 1], weights=[uncertainty_val, 1 - uncertainty_val], k=population_amount)

    # adding uncertainty to our model
    for k in range(population_amount):
        # person k has alleles j,l
        j, l = alleles_individuals[k]
        j, l = min(j, l), max(j, l)

        # choice whether this person will have uncertain alleles
        choice = choices[k]

        # this person has certain alleles
        if choice == 1:
            all_probabilities[j, l, k] = 1.0
        # this person has uncertain alleles
        if choice == 0:
            all_probabilities[:, :, k] = probabilities

    # now we have the probabilities {p_k(i,j)}
    # lets get the final p(i,j) = 1 / N * sum_k p_k(i,j)
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            probability = 0.0
            for k in range(population_amount):
                if t == m:
                    probability += all_probabilities[t, m, k]
                else:
                    probability += (all_probabilities[t, m, k] + all_probabilities[m, t, k])
            probability /= population_amount

            observed_probabilities[t, m] = probability

    observations = observed_probabilities * population_amount

    # we need to calculate alleles probabilities (from observations)
    alleles_probabilities = np.zeros(alleles_count)
    for i in range(alleles_count):
        for j in range(i, alleles_count):
            alleles_probabilities[i] += 0.5 * observed_probabilities[i, j]
            alleles_probabilities[j] += 0.5 * observed_probabilities[i, j]

    correction_vector = []
    variance_no_sampling_vector = []

    # correction = np.zeros(shape=(alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # instead of dividing by zero, we set the correction to 1
            if observations[j, k] == 0:
                correction = 1.0
            else:
                squared_sum = 0.0
                for l in range(population_amount):
                    squared_sum += (all_probabilities[j, k, l] ** 2)
                correction = squared_sum / observations[j, k]
            mult = 2
            if j == k:
                mult = 1
            expected = mult * population_amount * alleles_probabilities[j] * alleles_probabilities[k]
            correction_vector.append(expected * correction)

    observations_vector = []
    for t in range(alleles_count):
        for m in range(t, alleles_count):
            probabilities_list = all_probabilities[t, m, :].flatten()
            mean_value = np.mean(probabilities_list)
            variance = calc_variance(input_list=probabilities_list,
                                     input_mean=mean_value)
            variance_no_sampling_vector.append(variance)

            observations_vector.append(observations[t, m])
    list_values = variance_no_sampling_vector + observations_vector + correction_vector
    min_val = min(list_values)
    max_value = max(list_values)

    print(variance_no_sampling_vector)

    plt.scatter(variance_no_sampling_vector, observations_vector, label='OBSERVED', color='deeppink')
    plt.scatter(variance_no_sampling_vector, correction_vector, label='VAR CORRECTED', color='slategrey')
    plt.xlabel('VAR NO SAMPLING', fontsize=16)
    plt.xlim([min_val, max_value])
    plt.ylim([min_val, max_value])
    ax = plt.gca()
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    plt.legend(fontsize=12)


def prepare_probabilities(alleles_count):
    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)

    return alleles_probabilities


if __name__ == '__main__':
    alleles_amount = 50  # 20
    population_size = 10000  # 10000
    alpha_vals = [1.0]
    uncertainty_vals = [0.0, 0.2, 0.4]

    sns.set_style('white')

    plt.rcParams["font.family"] = "Arial"

    # fig, axes = plt.subplots(len(alpha_vals), len(uncertainty_vals), constrained_layout=True,
    #                          figsize=(10, 4))

    print('Calculating probabilities')
    # get {p(i)}, {p(i|j)}
    alleles_probabilities_ = prepare_probabilities(alleles_count=alleles_amount)

    plot_num = 0

    for alpha in alpha_vals:
        for uncertainty in uncertainty_vals:
            print(f'alpha: {alpha}, uncertainty: {uncertainty}')
            plot_num += 1
            plt.figure(figsize=((6, 6)))
            sns.set_style('white')
            plt.rcParams["font.family"] = "Arial"
            run_experiment(alleles_count=alleles_amount,
                           population_amount=population_size,
                           alpha_val=alpha, uncertainty_val=uncertainty,
                           alleles_probabilities=alleles_probabilities_,
                           plot_index=plot_num)
            uncertainty_str = str(int(uncertainty * 10))
            plt.savefig(f'variance_no_sampling_vs_observed_and_var_corrected_uncertainty_{uncertainty_str}.png', pad_inches=0.2, bbox_inches="tight")
    # plt.savefig('variance_no_sampling_vs_observed_and_var_corrected.png', pad_inches=0.2, bbox_inches="tight")
    # plt.show()