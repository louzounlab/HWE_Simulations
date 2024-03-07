from hwetests import asta
import utils

if __name__ == '__main__':
    # importing levels and 5 races names
    levels_set = utils.LEVELS_LIST
    races_set = utils.RACES_5_LIST

    for level in levels_set:
        for race in races_set:
            print(f'Current level: {level}, race: {race}')
            path = f'../../../data/hla_data_5_populations/levels/{level}/races/{race}/id_allele1_allele2_probability'
            file_to_save_name = f'plot_level_{level}_race_{race}.png'
            p_val, statistic, dof = asta.full_algorithm(file_path=path,
                                                        cutoff_value=2.0,
                                                        should_save_plot=file_to_save_name)
