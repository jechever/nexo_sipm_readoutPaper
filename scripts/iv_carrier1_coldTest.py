import matplotlib.pyplot as plt
import pandas as pd
#import numpy as np
"""
IV plotter for data taken by Gabriele at 165 K in probe station
"""

# cd Documents/nEXO/sipm_analysis/code/

# ----------------------------
# -------- Carrier 1 ---------
# ----------------------------

plt.rcParams.update({'font.size': 12})
plt.figure() 

data = pd.read_csv("../IV_data_carrierBoards/cold_test/ivs_165K.csv")

#data = pd.read_csv("../uiucArray_ivs/copy-s4_carrier1.csv")
current = abs(data[['ch1']])
volt = data[['Voltage2']]
plt.plot(abs(volt), current, label = 's1')

#data = pd.read_csv("../uiucArray_ivs/copy-s2_carrier1.csv")
current = abs(data[['ch2']])
volt = data[['Voltage1']]
plt.plot(abs(volt), current, label = 's2')

current = abs(data[['ch3']])
volt = data[['Voltage1']]
plt.plot(abs(volt), current, label = 's3')


#data = pd.read_csv("../uiucArray_ivs/copy-s3_carrier1.csv")
current = abs(data[['ch6']])
volt = data[['Voltage2']]
plt.plot(abs(volt), current, label = 's6')

plt.title('carrier 1')
plt.ylabel('$|I|$ (A)')
plt.xlabel('$|V_{bias}|$ (V)')
plt.legend()
plt.grid()
plt.show()

plt.yscale('log')

