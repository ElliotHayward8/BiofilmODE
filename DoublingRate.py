import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('CellsvsTime.csv')

data = pd.DataFrame(df).to_numpy()

i = 1
delta_N = []
while i < len(data) - 59:
    delta_N.append((data[i + 59][1]/data[i - 1][1]) - 1)
    i += 1

print(sum(delta_N))
print(len(delta_N))

print(sum(delta_N)/len(delta_N))


plt.plot(df['Minute'], df['Num cells'])
plt.grid()
plt.show()
