from models_simulated import chi_squared_second_attempt

file_path = 'simulation_data.csv'

data = []
with open(file_path, encoding="utf8") as infile:
    for index, line in enumerate(infile):
        # header
        if index == 0:
            continue
        lst = line.strip('\n').replace('+', ' ').replace(',', ' ').split()
        if lst[6] == '?':
            continue
        data.append([float(x) for x in lst])

chi_squared_second_attempt.plot_variance_vs_corrected_variance(data=data)
