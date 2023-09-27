############# CONSTANTS ############

# each race represents a different population
RACES_SET = {
    'CAU',
    'HIS',
    'AFA',
    'NAM',
    'API'
}
# there is also 'UNK', stands for unknown race.

# levels of the alleles
LEVELS_SET = {
    'A',
    'B',
    'C',
    'DQB1',
    'DRB1'
}

# dictionary of levels -> relative positions in the real data file
LEVELS_DICT = {'A': 0,
               'B': 1,
               'C': 2,
               'DQB1': 3,
               'DRB1': 4}

# dictionary of level -> race -> amount of people in race
AMOUNTS_DICT_DICTS = {'A': {'AFA': 653415,
                            'API': 804196,
                            'CAU': 4656668,
                            'HIS': 1125460,
                            'NAM': 57577},
                      'B': {'AFA': 653415,
                            'API': 804196,
                            'CAU': 4656668,
                            'HIS': 1125460,
                            'NAM': 57577},
                      'C': {'AFA': 653415,
                            'API': 804196,
                            'CAU': 4656668,
                            'HIS': 1125460,
                            'NAM': 57577},
                      'DQB1': {'AFA': 653415,
                               'API': 804196,
                               'CAU': 4656668,
                               'HIS': 1125460,
                               'NAM': 57577},
                      'DRB1': {'AFA': 653415,
                               'API': 804196,
                               'CAU': 4656668,
                               'HIS': 1125460,
                               'NAM': 57577}}

OLD_CHI_SQUARED = 'old_chi_squared'
NEW_CHI_SQUARED = 'new_chi_squared'

# 51,655,587 lines in us_grma
