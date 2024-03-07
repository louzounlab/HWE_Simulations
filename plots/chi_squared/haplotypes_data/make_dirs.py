import os
import constants

# THIS FILE CREATE THE DIRECTORIES USED FOR THE HAPLOTYPES DATA


def create_files(file):
    # path to this directory
    current_path = os.getcwd()

    # -> make directory: data -> real_data -> levels
    file_path = os.path.join(current_path, "../../../data", "haplotypes_data", file)
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    # races corresponding to this file
    races_list = constants.FILE_NAME_TO_POPULATIONS[file]
    # levels
    levels_list = constants.LEVELS_LIST
    # create file -> races
    races_path = os.path.join(file_path, "races")
    if not os.path.exists(races_path):
        os.makedirs(races_path)
    # creating races -> race
    for race in races_list:
        current_race_path = os.path.join(races_path, race)
        if not os.path.exists(current_race_path):
            os.makedirs(current_race_path)
    # create file -> levels
    levels_path = os.path.join(file_path, "levels")
    if not os.path.exists(levels_path):
        os.makedirs(levels_path)
    for level in levels_list:
        # for each level create directory: levels -> level
        current_level_path = os.path.join(levels_path, level)
        if not os.path.exists(current_level_path):
            os.makedirs(current_level_path)
        # -> races
        races_path = os.path.join(current_level_path, "races")
        # for each level and race create directory: levels -> level -> races -> race
        for race in races_list:
            current_race_path = os.path.join(races_path, race)
            if not os.path.exists(current_race_path):
                os.makedirs(current_race_path)
    # create file -> haplotypes
    haplotypes_path = os.path.join(file_path, "haplotypes")
    if not os.path.exists(haplotypes_path):
        os.makedirs(haplotypes_path)
    haplotypes_races_path = os.path.join(haplotypes_path, "races")
    if not os.path.exists(haplotypes_races_path):
        os.makedirs(haplotypes_races_path)
    # for each race in haplotypes create directory: haplotypes -> races -> race
    for race in races_list:
        current_race_path = os.path.join(haplotypes_races_path, race)
        if not os.path.exists(current_race_path):
            os.makedirs(current_race_path)
    # create results directory
    results_path = os.path.join(file_path, "results")
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    # create results -> haplotypes
    results_haplotypes_path = os.path.join(results_path, "haplotypes")
    if not os.path.exists(results_haplotypes_path):
        os.makedirs(results_haplotypes_path)
    # create results -> levels
    results_levels_path = os.path.join(results_path, "levels")
    if not os.path.exists(results_levels_path):
        os.makedirs(results_levels_path)


if __name__ == '__main__':
    files = constants.FILE_NAMES
    for file_ in files:
        create_files(file_)

