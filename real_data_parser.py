# [ID],[ALLELE1],[ALLELE2],[PROBABILITY],[AMOUNT]
# population: 8478533, alleles: 1487338

import numpy as np


# parse data into csv files. population_probs, id dictionary and alleles dictionary
def parse_data():
    pass


# get data from the csv files. returns: population_probs,
def get_data():
    pass

dict_id = {}
dict_alleles = {}

# alleles_count = 0
# population_amount = 0

with open("data/us_grma.haps.freqs", encoding="utf8") as infile:
    for index, line in enumerate(infile):
        lst = line.split(',')
        # print(line)
#         if len(lst) != 5:
#             # problem with the line
#             continue
#         # getting the information from the current line
#         id_str = lst[0]
#         allele_1 = lst[1]
#         allele_2 = lst[2]
#         probability = lst[3]
#         amount = lst[4]
#
#         # updating id dictionary for the id's
#         if id_str not in dict_id:
#             dict_id[id_str] = index
#             population_amount += 1
#
#         # updating dictionary for the alleles
#         for allele in (allele_1, allele_2):
#             if allele not in dict_alleles:
#                 dict_alleles[allele] = alleles_count
#                 alleles_count += 1

alleles_count = 1487338
population_amount = 8478533

population_probs = np.zeros(shape=(alleles_count, alleles_count, population_amount))

alleles_counter = 0
population_counter = 0

with open("data/us_grma.haps.freqs", encoding="utf8") as infile:
    for index, line in enumerate(infile):
        lst = line.split(',')
        if lst(lst) != 5:
            continue

        # getting information from current line
        id_str = lst[0]
        allele_1_str = lst[1]
        allele_2_str = lst[2]
        probability_str = lst[3]
        amount_str = lst[4]

        # updating id dictionary for the id's
        if id_str not in dict_id:
            dict_id[id_str] = population_counter
            population_counter += 1

        # updating dictionary for the alleles
        for allele in (allele_1, allele_2):
            if allele not in dict_alleles:
                dict_alleles[allele] = alleles_counter
                alleles_counter += 1

        # getting values of current line from dictionaries
        id = dict_id[id_str]
        allele_1 = dict_alleles[allele_1_str]
        allele_2 = dict_alleles[allele_2_str]
        amount = int(amount_str)

        population_probs[allele_1, allele_2, id] += amount
