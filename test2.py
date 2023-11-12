d = {'1_3': [1,2]}

for key, value in d.items():
    i_str, j_str = key.split('_')
    print(type(j_str))
    print(i_str, j_str)
    print(value)