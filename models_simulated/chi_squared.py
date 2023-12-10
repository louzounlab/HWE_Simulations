import random
import numpy as np
import scipy.stats as stats


#####
# Special comments:
# EVEN: means we take care of the case i!=j different from the case i=j
# DEBUG: didn't debug this part yet
#####


# assuming x is 2-dimensional
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


# data: N x 2, row = (j, k). alleles j,k
def count_row_occurrence_in_2d_array(j, k, data):
    count = 0
    # for every row (person) in population
    for row in range(data.shape[0]):
        if (data[row, 0] == j and data[row, 1] == k) or (data[row, 0] == k and data[row, 1] == j):
            count += 1
            # count += (alleles_probabilities[j] * alleles_probabilities[k] / (population_probs[j,k,row]))
            # count += 1

    return count


# need to calculate only on the upper triangle because the matrices are symmetric
# test_type: 0 - old, 1 - correction, 2 - sampling
def calculate_chi_squared_value(counts_expected_, counts_observed_, variances, correction_, test_type):
    value = 0
    for row in range(counts_expected_.shape[0]):
        for col in range(row, counts_expected_.shape[1]):
            expected_ = counts_expected_[row, col]
            observed_ = counts_observed_[row, col]
            variance_ = variances[row, col]
            if test_type in {0, 2}:
                value += (((expected_ - observed_) ** 2) / variance_)
            else:
                # value += correction_[row, col] * (((expected_ - observed_) ** 2) / variance_)
                value += (1 / correction_[row, col]) * (((expected_ - observed_) ** 2) / variance_)
    return value


def calc_variance(input_list, input_mean):
    var = 0.0
    for element in input_list:
        var += (element - input_mean) ** 2
    return var


def prepare_probabilities(alleles_count):
    # probabilities {p(i)}
    alleles_probabilities = calculate_alleles_probabilities(alleles_count)
    return alleles_probabilities


def run_experiment(alleles_count, population_amount, alpha_val, uncertainty_val,
                   alleles_probabilities):
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

    # now we multiply the upper triangle by 2
    # for t in range(alleles_count):
    #     for m in range(t + 1, alleles_count):
    #         probabilities[t, m] *= 2

    # print(f'alleles_individuals: {alleles_individuals}')

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

    # sum_observed = 0.0
    # for k in range(alleles_count):
    #     for j in range(k, alleles_count):
    #         sum_observed += observations[k, j]

    counts_expected = np.zeros((alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # EVEN
            mult = 1
            if k != j:
                mult = 2
            expected_value = mult * population_amount * alleles_probabilities[k] * alleles_probabilities[j]
            counts_expected[k, j] = expected_value
            counts_expected[j, k] = expected_value

    correction = np.zeros(shape=(alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # instead of dividing by zero, we set the correction to 1
            if observations[j, k] == 0:
                correction[j, k] = 1.0
            else:
                squared_sum = 0.0
                for l in range(population_amount):
                    squared_sum += (all_probabilities[j, k, l] ** 2)
                correction[j, k] = squared_sum / observations[j, k]

    variances_ = np.zeros(shape=(alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # EVEN
            mult = 1
            if k != j:
                mult = 2
            expected_value = mult * population_amount * alleles_probabilities[k] * alleles_probabilities[j]
            variances_[j, k] = expected_value

    chi_squared_stat_old = calculate_chi_squared_value(counts_expected_=counts_expected,
                                                       counts_observed_=observations,
                                                       variances=variances_,
                                                       correction_=correction,
                                                       test_type=0)

    chi_squared_stat_corrected = calculate_chi_squared_value(counts_expected_=counts_expected,
                                                             counts_observed_=observations,
                                                             variances=variances_,
                                                             correction_=correction,
                                                             test_type=1)

    # for each k make {p_k(i,j)} symmetric and choose for every k alleles i,j
    for k in range(population_amount):
        # for t in range(alleles_count):
        #     for m in range(t, alleles_count):
        #         all_probabilities[t, m, k] = (all_probabilities[t, m, k] + all_probabilities[m, t, k]) / 2
        probabilities_list = all_probabilities[:, :, k].flatten()
        index = random.choices(population=range(len(probabilities_list)), weights=probabilities_list, k=1)[0]

        # the right allele of this index element
        col = index % alleles_count
        # the left allele
        row = (index - col) // alleles_count
        # making sure i <= j
        row, col = min(row, col), max(row, col)
        # assignment of the alleles to the person
        alleles_individuals[k, :] = [row, col]

    for t in range(alleles_count):
        for m in range(t, alleles_count):
            probabilities_list = all_probabilities[t, m, :].flatten()
            mean_value = np.mean(probabilities_list)
            variance = calc_variance(input_list=probabilities_list,
                                     input_mean=mean_value)
            variances_[t, m] = variance
    # print(variances_)

    chi_squared_stat_sampling = calculate_chi_squared_value(counts_expected_=counts_expected,
                                                            counts_observed_=observations,
                                                            variances=variances_,
                                                            correction_=correction,
                                                            test_type=2)

    dof = (alleles_count * (alleles_count + 1)) / 2 - 1

    # print(f' alpha for choice: {alpha_val}')
    # print(f' chi square value: {chi_squared_stat}')

    # crit = stats.chi2.ppf(q=0.95, df=dof)
    # print(f'Critical value: {crit}')

    p_value_old = 1 - stats.chi2.cdf(x=chi_squared_stat_old,
                                     df=dof)
    p_value_corrected = 1 - stats.chi2.cdf(x=chi_squared_stat_corrected,
                                           df=dof)
    p_value_sampling = 1 - stats.chi2.cdf(x=chi_squared_stat_sampling,
                                          df=dof)
    return int(p_value_old < 0.05), int(p_value_corrected < 0.05), int(p_value_sampling < 0.05)
