import matplotlib.pyplot as plt
import pandas as pd
#import numpy as np
"""
IVs for data at RT originally taken by Gerard right after board assembly
"""

# cd Documents/nEXO/sipm_analysis/code/

# ----------------------------
# -------- Carrier 1 ---------
# ----------------------------

plt.rcParams.update({'font.size': 12})
plt.figure() 

data = pd.read_csv("../IV_data/carrier1_sipm1_run2.log")
i_array = data[['current']]
volt = data[['voltage']]
array = []
for index, row in volt.iterrows():
    val = row['voltage']
    new_val = val[0:8]
    array.append(abs(float(new_val)))
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[4:16]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's1')
plt.title('carrier 1')


data = pd.read_csv("../IV_data/carrier1_sipm2_run2.log")
i_array = data[['current']]
volt = data[['voltage']]
array = []
for index, row in volt.iterrows():
    val = row['voltage']
    new_val = val[0:8]
    array.append(abs(float(new_val)))
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[4:16]
    array.append(new_val)
    #print(new_val)
#array1 = array[::-1]
plt.plot(x, array, label = 's2')


data = pd.read_csv("../IV_data/carrier1_sipm3_run2.log")
i_array = data[['current']]
volt = data[['voltage']]
array = []
for index, row in volt.iterrows():
    val = row['voltage']
    new_val = val[0:8]
    array.append(abs(float(new_val)))
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[4:16]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's3')


data = pd.read_csv("../IV_data/carrier1_sipm4_run2.log")
i_array = data[['current']]
volt = data[['voltage']]
array = []
for index, row in volt.iterrows():
    val = row['voltage']
    new_val = val[0:8]
    array.append(abs(float(new_val)))
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[4:16]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's4')


data = pd.read_csv("../IV_data/carrier1_sipm5_run2.log")
i_array = data[['current']]
volt = data[['voltage']]
array = []
for index, row in volt.iterrows():
    val = row['voltage']
    new_val = val[0:8]
    array.append(abs(float(new_val)))
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[4:16]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's5')


data = pd.read_csv("../IV_data/carrier1_sipm6_run2.log")
i_array = data[['current']]
volt = data[['voltage']]
array = []
for index, row in volt.iterrows():
    val = row['voltage']
    new_val = val[0:8]
    array.append(abs(float(new_val)))
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[4:16]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's6')

plt.ylabel('$|I|$ (A)')
plt.xlabel('$|V_{bias}|$ (V)')
plt.legend()
plt.grid()
plt.show()

plt.yscale('log')

# ----------------------------
# -------- Carrier 2 ---------
# ----------------------------

plt.figure() 
plt.rcParams.update({'font.size': 12})

data = pd.read_csv("../IV_data/carrier2_sipm1_run3.log")
i_array = data[['current']]
#x = np.linspace(0,50,len(i_array))
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[3:15]
    array.append(new_val)
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[20:32]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's1')
plt.title('carrier 2')


data = pd.read_csv("../IV_data/carrier2_sipm3_run2.log")
i_array = data[['current']]
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[3:15]
    array.append(new_val)
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[20:32]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's2')


data = pd.read_csv("../IV_data/carrier2_sipm3_run2.log")
i_array = data[['current']]
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[3:15]
    array.append(new_val)
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[20:32]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's3')


data = pd.read_csv("../IV_data/carrier2_sipm4_run2.log")
i_array = data[['current']]
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[3:15]
    array.append(new_val)
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[20:32]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's4')


data = pd.read_csv("../IV_data/carrier2_sipm5_run2.log")
i_array = data[['current']]
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[3:15]
    array.append(new_val)
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[20:32]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's5')


data = pd.read_csv("../IV_data/carrier2_sipm6_run2.log")
i_array = data[['current']]
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[3:15]
    array.append(new_val)
    #print(new_val)
x = array 
array = []
for index, row in i_array.iterrows():
    val = row['current']
    new_val = val[20:32]
    array.append(new_val)
    #print(new_val)
plt.plot(x, array, label = 's6')

plt.ylabel('$|I|$ (A)')
plt.xlabel('$|V_{bias}|$ (V)')
plt.legend()
plt.grid()
plt.show()

plt.yscale('log')