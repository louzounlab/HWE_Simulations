############# CONSTANTS ############

# each race represents a different population
RACES_5_LIST = [
    'CAU',
    'HIS',
    'AFA',
    'NAM',
    'API'
]

RACES_21_LIST = ["AAFA",
                 "AFB",
                 "AINDI",
                 "AISC",
                 "ALANAM",
                 "AMIND",
                 "CARB",
                 "CARHIS",
                 "CARIBI",
                 "HAWI",
                 "FILII",
                 "KORI",
                 "JAPI",
                 "MSWHIS",
                 "MENAFC",
                 "NAMER",
                 "NCHI",
                 "SCAHIS",
                 "SCAMB",
                 "SCSEAI",
                 "VIET"]

RACES_1_LIST = ['ALL']

# there is also 'UNK', stands for unknown race.

# levels of the alleles
LEVELS_LIST = [
    'A',
    'B',
    'C',
    'DQB1',
    'DRB1'
]

# dictionary: small populations -> broad populations
races_to_broad_races = {"AAFA": "AFA",
                        "AFB": "AFA",
                        "AINDI": "API",
                        "AISC": "NAM",
                        "ALANAM": "NAM",
                        "AMIND": "NAM",
                        "CARB": "AFA",
                        "CARHIS": "HIS",
                        "CARIBI": "NAM",
                        "HAWI": "API",
                        "FILII": "API",
                        "KORI": "API",
                        "JAPI": "API",
                        "MSWHIS": "HIS",
                        "MENAFC": "CAU",
                        "NAMER": "CAU",
                        "NCHI": "API",
                        "SCAHIS": "HIS",
                        "SCAMB": "AFA",
                        "SCSEAI": "API",
                        "VIET": "API"}
broad_races_to_colors = {'CAU': 'magenta',
                         'HIS': 'b',
                         'AFA': 'r',
                         'NAM': 'g',
                         'API': 'k'}

# dictionary of levels -> relative positions in the real data file
LEVELS_DICT = {'A': 0,
               'B': 1,
               'C': 2,
               'DQB1': 3,
               'DRB1': 4}

# # dictionary of level -> race -> amount of people in race
# AMOUNTS_5_DICT_DICTS = {'A': {'AFA': 653415,
#                               'API': 804196,
#                               'CAU': 4656668,
#                               'HIS': 1125460,
#                               'NAM': 57577},
#                         'B': {'AFA': 653415,
#                               'API': 804196,
#                               'CAU': 4656668,
#                               'HIS': 1125460,
#                               'NAM': 57577},
#                         'C': {'AFA': 653415,
#                               'API': 804196,
#                               'CAU': 4656668,
#                               'HIS': 1125460,
#                               'NAM': 57577},
#                         'DQB1': {'AFA': 653415,
#                                  'API': 804196,
#                                  'CAU': 4656668,
#                                  'HIS': 1125460,
#                                  'NAM': 57577},
#                         'DRB1': {'AFA': 653415,
#                                  'API': 804196,
#                                  'CAU': 4656668,
#                                  'HIS': 1125460,
#                                  'NAM': 57577}}
#
# AMOUNTS_21_DICT_DICTS = {
#     'VIET': 63958,
#     'KORI': 94604,
#     'FILII': 71565,
#     'CARB': 46793,
#     'HAWI': 15242,
#     'ALANAM': 4052,
#     'AINDI': 265827,
#     'CARIBI': 3113,
#     'SCAMB': 5430,
#     'NAMER': 2407199,
#     'SCSEAI': 46370,
#     'SCAHIS': 201472,
#     'MSWHIS': 77323,
#     'MENAFC': 130269,
#     'AAFA': 493089,
#     'AMIND': 35585,
#     'AISC': 8153,
#     'NCHI': 143571,
#     'AFB': 44218,
#     'CARHIS': 131405,
#     'JAPI': 22894
# }

OLD_CHI_SQUARED = 'old_chi_squared'
NEW_CHI_SQUARED = 'new_chi_squared'

# 51,655,587 lines in us_grma
