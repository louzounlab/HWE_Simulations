import random
import numpy as np
import scipy.stats as stats


alleles_count = 0
max_allele = 0

cut_off_value = 4.0
file_path = 'simulation_data.csv'
# first read all the rows and get indices of ids and alleles and amounts
with open(file_path, encoding="utf8") as infile:
    for index, line in enumerate(infile):
        # header
        if index == 0:
            continue
        lst = line.strip('\n').replace('+', ' ').replace(',', ' ').split()
        i = int(lst[0])
        j = int(lst[1])
        max_allele = max(max_allele, i, j)
alleles_count = max_allele + 1

observed = np.zeros((alleles_count, alleles_count))
expected = np.zeros((alleles_count, alleles_count))
rho = np.zeros((alleles_count, alleles_count))

with open(file_path, encoding="utf8") as infile:
    for index, line in enumerate(infile):
        # header
        if index == 0:
            continue
        lst = line.strip('\n').replace('+', ' ').replace(',', ' ').split()
        i = int(lst[0])
        j = int(lst[1])
        i, j = min(i, j), max(i, j)
        observed[i, j] = float(lst[2])
        expected[i, j] = float(lst[3])
        rho[i, j] = float(lst[4])
statistic = 0.0
amount_of_small_expected = 0
couples_amount = (alleles_count * (alleles_count + 1)) / 2 - 1
print(observed)
for i in range(alleles_count):
    for j in range(i, alleles_count):
        if expected[i, j] < cut_off_value:
            amount_of_small_expected += 1
            continue
        statistic += ((observed[i, j] - expected[i, j]) ** 2 / ((expected[i, j]) * rho[i, j]))

dof_old = couples_amount - amount_of_small_expected
print(dof_old)

p_value_old = 1 - stats.chi2.cdf(x=statistic,
                                 df=dof_old)
print(p_value_old)
print(statistic)
