from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np 

'''
Script for producing reso plots as function of sipm N for 2019 sipm paper 
'''

# -----------------------------------------------------------------------------
# -----------------------------   SERIES --------------------------------------
# -----------------------------------------------------------------------------

fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
# Series 

resolutions =[0.063, 0.134, 0.17] ; nsipms = [2, 4, 6]; bias = '3.5';
m,b = np.polyfit(nsipms,resolutions,1)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit)
#plt.figure(); 
plt.plot(nsipms, resolutions,'co'); plt.grid(); plt.plot(nsipms, yList, 'c', linestyle='dashed', label ='$V_{over}$ = '+bias+'V');

resolutions =[0.072, 0.124, 0.17] ; nsipms = [2, 4, 6]; bias = '3.5';
m,b = np.polyfit(nsipms,resolutions,1)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit)
#plt.figure(); 
plt.plot(nsipms, resolutions,'mo'); plt.grid(); plt.plot(nsipms, yList, 'm', linestyle='dashed', label ='$V_{over}$ = '+bias+'V');

resolutions =[0.0781, 0.101, 0.129] ; nsipms = [2, 4, 6]; bias = '4.5';
m,b = np.polyfit(nsipms,resolutions,1)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit)
plt.plot(nsipms, resolutions,'ko'); plt.plot(nsipms, yList, 'k',linestyle='dashed', label ='$V_{over}$ = '+bias+'V');
resolutions =[0.0738, 0.096, 0.125] ; nsipms = [2, 4, 6]; bias = '5.0';
m,b = np.polyfit(nsipms,resolutions,1)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit) 
plt.plot(nsipms, resolutions,'go'); plt.plot(nsipms, yList, 'g', linestyle='dashed', label ='$V_{over}$ = '+bias+'V'); 
#plt.axhline(y=0.1, color='r', linestyle='-', label = 'nEXO requirement')
plt.xlabel('Number of SiPMs'); plt.ylabel('$\sigma_{SPE}$'); plt.show(); plt.legend();
plt.grid(); plt.xlim([1.5, 6.2])
plt.ylim([0, 0.25])

resolutions =[0.079, 0.081, 0.1] ; nsipms = [2, 4, 6]; bias = '5.5';
m,b = np.polyfit(nsipms,resolutions,1)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit) 
plt.plot(nsipms, resolutions,'bo'); plt.plot(nsipms, yList, 'b', linestyle='dashed', label ='$V_{over}$ = '+bias+'V'); 
#plt.axhline(y=0.1, color='r', linestyle='-', label = 'nEXO requirement')
plt.xlabel('Number of SiPMs'); plt.ylabel('$\sigma_{SPE}$'); 
plt.grid(); plt.show(); plt.legend(); 
plt.xlim([1.8, 6.2]); plt.ylim([0.05, 0.25])

# -----------------------------------------------------------------------------
# -----------------------------   PARALLEL ------------------------------------
# -----------------------------------------------------------------------------

resolutions =[0.118, 0.164, 0.201, 0.236] ; nsipms = [3, 4, 5, 6]; bias = '4.5';
m,b = np.polyfit(nsipms,resolutions,1)
fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
print('m = ',m)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit) 
#plt.figure()
plt.plot(nsipms, resolutions,'ko'); plt.plot(nsipms, yList, 'k', linestyle='dashed', label ='$V_{over}$ = '+bias+'V'); 

resolutions =[0.103, 0.157, 0.17, 0.198] ; nsipms = [3, 4, 5, 6]; bias = '5.0';
m,b = np.polyfit(nsipms,resolutions,1)
print('m = ',m)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit)
plt.plot(nsipms, resolutions,'go'); plt.plot(nsipms, yList, 'g',linestyle='dashed', label ='$V_{over}$ = '+bias+'V');

resolutions =[0.099, 0.144, 0.165, 0.18] ; nsipms = [3, 4, 5, 6]; bias = '5.5';
m,b = np.polyfit(nsipms,resolutions,1)
print('m = ',m)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit)
plt.plot(nsipms, resolutions,'bo'); plt.plot(nsipms, yList, 'b', linestyle='dashed', label ='$V_{over}$ = '+bias+'V');
plt.xlabel('Number of SiPMs'); plt.ylabel('$\sigma_{SPE}$'); 
plt.grid(); plt.show(); plt.legend(); 
plt.xlim([1.8, 6.2]); plt.ylim([0.05, 0.25])