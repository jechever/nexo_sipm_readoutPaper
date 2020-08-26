import matplotlib.pyplot as plt
import pandas as pd
#import numpy as np
"""
IV plotter for data taken by Gabriele at RT (295 K) in probe station
"""

# cd Documents/nEXO/sipm_analysis/code/

# ----------------------------
# -------- Carrier 1 ---------
# ----------------------------

plt.rcParams.update({'font.size': 12})
plt.figure() 

data = pd.read_csv("../IV_data_carrierBoards/rt_gabriele/sipm1_adj.csv")
current = abs(data[[' IK']])
volt = data[[' VA']]
plt.plot(abs(volt), current, label = 's1')

data = pd.read_csv("../IV_data_carrierBoards/rt_gabriele/sipm2_adj.csv")
current = abs(data[[' IK']])
volt = data[[' VA']]
plt.plot(abs(volt), current, label = 's2')

data = pd.read_csv("../IV_data_carrierBoards/rt_gabriele/sipm3_adj.csv")
current = abs(data[[' IK']])
volt = data[[' VA']]
plt.plot(abs(volt), current, label = 's3')

data = pd.read_csv("../IV_data_carrierBoards/rt_gabriele/sipm4_adj.csv")
current = abs(data[[' IK']])
volt = data[[' VA']]
plt.plot(abs(volt), current, label = 's4')

data = pd.read_csv("../IV_data_carrierBoards/rt_gabriele/sipm5_adj.csv")
current = abs(data[[' IK']])
volt = data[[' VA']]
plt.plot(abs(volt), current, label = 's5')

data = pd.read_csv("../IV_data_carrierBoards/rt_gabriele/sipm6_adj.csv")
current = abs(data[[' IK']])
volt = data[[' VA']]
plt.plot(abs(volt), current, label = 's6')

# Edit axes properties 
plt.title('carrier 1')
plt.ylabel('$|I|$ (A)')
plt.xlabel('$|V_{bias}|$ (V)')
plt.legend()
plt.grid()
plt.show()

plt.yscale('log')

