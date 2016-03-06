from scipy.io import loadmat

from mutate import *


hiv_data = loadmat('flhivdata.mat')


# print(hiv_data.keys())

# print(hiv_data["lc5"])

# print(hiv_data["lc1"])

# print(hiv_data["__globals__"])


print('tree sim')
print(evolution_simulator(.1, 100000, 10, 'GATTACA', simple_killing, 2))






