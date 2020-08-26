"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 4-SiPMs 
(2s2p) for t-rise dependence study, using Gaussian filter. Uses sample biased 
at Vbias = -68 V (Vover = 4.5 V)
"""
from sipm_signalProcessing import getAvgPulse_gauss, get_amplitudes_noAP, plot_filteredPulses_gauss
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import detect_peaks
import fits
 
# SET PATH: 
# cd Documents/nEXO/sipm_analysis/code

# Read data file from Mietek: 
#df_parallel_gauss = pd.read_csv('simulations_mietek/enc_semigauss_shaper_parallel.csv')
df_series_gauss = pd.read_csv('simulations_mietek/enc_semigauss_shaper_2_in_series.csv')
tpeak_series_gauss = df_series_gauss['tpeak'].tolist()
enc_series_gauss_4sipm = df_series_gauss['4sipm_enc'].tolist()

# Initiate lists for observables:
resolutions_4s, Tlist_4s = [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; peak = 'single';
peMin_unfiltered = 0.02; 

# ------------------------------------------------------------------------------
# ----------------------------- 4-sipms 2s2p  ----------------------------------
# ------------------------------------------------------------------------------
sipmN = '4'; source = 'dark';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
"""
# --------------------------------- 50 ns  -------------------------------------

# Set peaking time:  
peMin = 0.017; peMax = 0.032; peakingTime = 50

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
    
# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 28; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));                  
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)
"""
## --------------------------------- 80 ns  -------------------------------------

# Set peaking time:  
peMin = 0.017; peMax = 0.028; peakingTime = 80

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 15, 25; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 100 ns  -------------------------------------

# Set peaking time:  
peMin = 0.015; peMax = 0.027; peakingTime = 100

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 12, 23; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 110 ns  -------------------------------------

# Set peaking time:  
peMin = 0.012; peMax = 0.025; peakingTime = 120

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
    
# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 12, 22; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 130 ns  -------------------------------------

# Set peaking time:  
peMin = 0.01; peMax = 0.024; peakingTime = 140

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 11, 20; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 140 ns  -------------------------------------

# Set peaking time:  
peMin = 0.01; peMax = 0.022; peakingTime = 160

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 10, 19; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 170 ns  -------------------------------------

# Set peaking time:  
peMin = 0.01; peMax = 0.02; peakingTime = 180

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 9, 18; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 180 ns  -------------------------------------

# Set peaking time:  
peMin = 0.01; peMax = 0.019; peakingTime = 200

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 9, 17; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 220 ns  -------------------------------------

# Set peaking time:  
peMin = 0.009; peMax = 0.018; peakingTime = 240

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 80; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
#plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 12, 22; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 240 ns  -------------------------------------

# Set peaking time:  
peMin = 0.009; peMax = 0.016; peakingTime = 280

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 60; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[1]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 14, 27; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 280 ns  -------------------------------------

# Set peaking time:  
peMin = 0.007; peMax = 0.015; peakingTime = 320

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 60; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[1]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 14, 25; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 330 ns  -------------------------------------

# Set peaking time:  
peMin = 0.006; peMax = 0.013; peakingTime = 400

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 40; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[1]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 32; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 400 ns  -------------------------------------

# Set peaking time:  
peMin = 0.005; peMax = 0.01; peakingTime = 520

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 40; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[1]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 13, 27; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; #delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 550 ns  -------------------------------------

# Set peaking time:  
peMin = 0.003; peMax = 0.009; peakingTime = 680

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 40; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 8, 22; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 620 ns  -------------------------------------

# Set peaking time:  
peMin = 0.003; peMax = 0.008; peakingTime = 780

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()


# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 10, 26; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 680 ns  -------------------------------------

# Set peaking time:  
peMin = 0.002; peMax = 0.007; peakingTime = 900

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()


# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 9, 23; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 870 ns  -------------------------------------

# Set peaking time:  
peMin = 0.001; peMax = 0.006; peakingTime = 1100

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 25; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
x = x/1.14;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 7, 24; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 870 ns  -------------------------------------

# Set peaking time:  
peMin = 0.001; peMax = 0.005; peakingTime = 1300

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()


# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 2, 19; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 1100 ns  -------------------------------------

# Set peaking time:  
peMin = 0.001; peMax = 0.005; peakingTime = 1600

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 20; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
x = x/1.1

# EXPECTED VALUES 
# Single PE
mu1 = 1; 
A1 = 0.3; 
# Second PE
mu2 = 2; 
A2 = 0.02;  
sig = 0.1; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()


# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 2, 21; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
plt.close() 

xAxis = np.linspace(0,2,50)
gaussianFit = fits.gauss(x, mu1, sig1, A1)
gaussianFit = list(gaussianFit)
index = gaussianFit.index(max(gaussianFit))
x = x/(x[index])
gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)

fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; 
print('RESOLUTION: '+str(sig1))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
#plt.figure(); plt.plot(y)

plt.close('all')

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs Tpeak plot: 
plt.figure(); 
Tlist_4s = [x*1e-9 for x in Tlist_4s] 
plt.plot(Tlist_4s, resolutions_4s, 'g.', label = 'Data')
plt.plot(tpeak_series_gauss, enc_series_gauss_4sipm, 'g', label = 'Simulations')
plt.title('Gauss filter, 4 SiPMs: 2s2p, $V_{over}$ = 4.5 V') ; 
plt.xlabel('$t_{peak}$ (s)'); plt.ylabel('enc (s.p.e.)');
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()
plt.savefig('enc-gaussFilter-noAP-4s.png')

# Make pandas dataframes: 
df_freq_4s = pd.DataFrame({'peaking time (s)': Tlist_4s,
                   'resolution (s.p.e.)': resolutions_4s})
                   
# Save as .csv file: 
df_freq_4s.to_csv('enc-gaussFilter-noAP-4s.csv')