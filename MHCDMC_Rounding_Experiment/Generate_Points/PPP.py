import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from random import randint
import scipy

# Simulation window parameters

x1Min = 0
x1Max = 800
y1Min = 0
y1Max = 800
x1Delta = x1Max - x1Min
y1Delta = y1Max - y1Min  # rectangle dimensions

x2Min = 1200
x2Max = 2000
y2Min = 0
y2Max = 800
x2Delta = x2Max - x2Min
y2Delta = y2Max - y2Min  # rectangle dimensions

x3Min = 600
x3Max = 1400
y3Min = 1200
y3Max = 2000
x3Delta = x3Max - x3Min
y3Delta = y3Max - y3Min  # rectangle dimensions

areaTotal = x1Delta * y1Delta

# Point process parameters
lambda0 = 200e-6;  # intensity (ie mean density) of the Poisson process
lambda1 = 5e-6
lambda2 = 50e-6

# Simulate Poisson point process
ms_num = scipy.stats.poisson(lambda0 * areaTotal).rvs()  # Poisson number of points
num_macro = scipy.stats.poisson(lambda1 * areaTotal).rvs()  # Poisson number of points
num_small = scipy.stats.poisson(lambda2 * areaTotal).rvs()  # Poisson number of points

print("ms_num : ", ms_num)
print("num_macro : ", num_macro)
print("num_small : ", num_small)

# coordinates of Poisson points
macro_X = x1Delta * scipy.stats.uniform.rvs(0, 1, (num_macro, 1)) + x1Min
macro_Y = y1Delta * scipy.stats.uniform.rvs(0, 1, (num_macro, 1)) + y1Min

small_X = x1Delta * scipy.stats.uniform.rvs(0, 1, (num_small, 1)) + x1Min
small_Y = y1Delta * scipy.stats.uniform.rvs(0, 1, (num_small, 1)) + y1Min

ms_X = x1Delta * scipy.stats.uniform.rvs(0, 1, (ms_num, 1)) + x1Min
ms_Y = y1Delta * scipy.stats.uniform.rvs(0, 1, (ms_num, 1)) + y1Min

# figure
plt.scatter(ms_X, ms_Y, s=1, c="black", marker="s", label='UE')
plt.xlabel("x");
plt.ylabel("y")

plt.scatter(macro_X, macro_Y, c="red", marker="^", label='Macro Cell')
plt.xlabel("x");
plt.ylabel("y")

plt.scatter(small_X, small_Y, c="green", marker="^", label='Small Cell')
plt.xlabel("x");
plt.ylabel("y")

plt.legend(bbox_to_anchor=(1, 1))
plt.show()