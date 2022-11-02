import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

df = pd.read_csv('CellsvsTime.csv')

data = pd.DataFrame(df).to_numpy()
data = data.astype(float)

num_cells = df['Num cells'].to_numpy()

dy = np.diff(num_cells)

idx = np.argwhere(dy)

data_fixed = num_cells.copy()
data_fixed = data_fixed.astype(float)

x = np.arange(num_cells.shape[0])


data_fixed[idx + 1] += 0.001

f = interp1d(data_fixed, x, fill_value='extrapolate')

data_half = num_cells/2.0

x_interp = f(data_half)

dbl = x - x_interp
# print(dbl)

hidden_dbl = np.ma.masked_invalid(dbl)

print(hidden_dbl)

print(hidden_dbl.sum()/(len(hidden_dbl) - 31))

# i = 1
# delta_N = []
# while i < len(data) - 59:
#     delta_N.append((data[i + 59][1]/data[i - 1][1]) - 1)
#     i += 1
#
# print(delta_N)
# print(len(delta_N))
#
# print(sum(delta_N)/len(delta_N))
#
#
# plt.plot(df['Minute'], df['Num cells'])
# plt.grid()
# plt.show()
