import os
import utils

# path to this directory
current_path = os.getcwd()

# -> make directory: data -> real_data -> levels
levels_path = os.path.join(current_path, "data", "hla_data_21_populations", "levels")
if not os.path.exists(levels_path):
    os.makedirs(levels_path)

# get levels and races from utils file
levels_set = utils.LEVELS_SET
races_set = utils.RACES_21_SET

# create data -> real_data_ races -> {race}
races_path = os.path.join(current_path, "data", "hla_data_21_populations", "races")
if not os.path.exists(races_path):
    os.makedirs(races_path)
for race in races_set:
    current_race_path = os.path.join(races_path, race)
    if not os.path.exists(current_race_path):
        os.makedirs(current_race_path)

for level in levels_set:
    # for each level create directory: levels -> level
    current_level_path = os.path.join(levels_path, level)
    if not os.path.exists(current_level_path):
        os.makedirs(current_level_path)

    # -> races
    races_path = os.path.join(current_level_path, "races")
    # for each level and race create directory: levels -> level -> races -> race
    for race in races_set:
        current_race_path = os.path.join(races_path, race)
        if not os.path.exists(current_race_path):
            os.makedirs(current_race_path)

# make directory for real data plots
plots_path = os.path.join(current_path, "plots")

chi_squared_tests = {utils.OLD_CHI_SQUARED, utils.NEW_CHI_SQUARED}

for test in chi_squared_tests:
    test_path = os.path.join(plots_path, test)
    if not os.path.exists(test_path):
        os.makedirs(test_path)


