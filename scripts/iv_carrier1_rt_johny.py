import matplotlib.pyplot as plt
import pandas as pd
#import numpy as np


# cd Documents/nEXO/sipm_analysis/code/

# ----------------------------
# -------- Carrier 1 ---------
# ----------------------------

plt.rcParams.update({'font.size': 12})
plt.figure() 

data = pd.read_csv("../uiucArray_ivs/copy-s1_carrier1.csv")
current = data[[' Io']]
volt = data[[' Vi']]
plt.plot(abs(volt), current, label = 's1')
plt.title('carrier 1')

data = pd.read_csv("../uiucArray_ivs/copy-s2_carrier1.csv")
current = data[[' Io']]
volt = data[[' Vi']]
plt.plot(abs(volt), current, label = 's2')

data = pd.read_csv("../uiucArray_ivs/copy-s3_carrier1.csv")
current = data[[' Io']]
volt = data[[' Vi']]
plt.plot(abs(volt), current, label = 's3')

data = pd.read_csv("../uiucArray_ivs/copy-s4_carrier1.csv")
current = data[[' Io']]
volt = data[[' Vi']]
plt.plot(abs(volt), current, label = 's4')

data = pd.read_csv("../uiucArray_ivs/copy-s5_carrier1.csv")
current = data[[' Io']]
volt = data[[' Vi']]
plt.plot(abs(volt), current, label = 's5')

data = pd.read_csv("../uiucArray_ivs/copy-s6_carrier1.csv")
current = data[[' Io']]
volt = data[[' Vi']]
plt.plot(abs(volt), current, label = 's6')

plt.ylabel('$|I|$ (A)')
plt.xlabel('$|V_{bias}|$ (V)')
plt.legend()
plt.grid()
plt.show()

plt.yscale('log')

