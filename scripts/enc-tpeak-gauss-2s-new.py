"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 2-SiPMs 
(2s) for t-rise dependence study, using Gaussian filter. Uses sample biased 
at Vbias = -68 V (Vover = 4.5 V)
"""
from sipm_signalProcessing import getAvgPulse_gauss, getAmplitudes_gauss, plot_filteredPulses_gauss
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import detect_peaks
import fits
 
# SET PATH: 
# cd Documents/nEXO/sipm_analysis/code

# Read data file from Mietek: 
df_series_gauss = pd.read_csv('simulations_mietek/enc_semigauss_shaper_2_in_series.csv')
tpeak_series_gauss = df_series_gauss['tpeak'].tolist()
enc_series_gauss_2sipm = df_series_gauss['2sipm_enc'].tolist()

# Initiate lists for observables:
resolutions_2s, Tlist_2s, errors = [], [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; peak = 'single';
peMin_unfiltered = 0.025;

# ------------------------------------------------------------------------------
# ----------------------------- 2-sipms 1s1p  ----------------------------------
# ------------------------------------------------------------------------------
sipmN = '2'; source = 'dark';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# --------------------------------- 50 ns  -------------------------------------
"""
# Set peaking time:  
peMin = 0.028; peMax = 0.043; peakingTime = 50

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True); 
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

# Plot PE spectrum:
lower, upper = 0, 160; nbins = 100; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 19, 27; 
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

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)
"""
# --------------------------------- 80 ns  -------------------------------------

# Set peaking time:  
peMin = 0.026; peMax = 0.04; peakingTime = 80

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

# Plot PE spectrum:
lower, upper = 0, 160; nbins = 100; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 17, 25; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 100 ns  -------------------------------------

# Set peaking time:  
peMin = 0.02; peMax = 0.038; peakingTime = 100

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); #plt.close()

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 10, mph = 0.02) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 21, 31; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 120 ns  -------------------------------------

# Set peaking time:  
peMin = 0.02; peMax = 0.036; peakingTime = 120

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); #plt.close()

# Plot PE spectrum:
lower, upper = 0, 140; nbins = 100; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 25; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 140 ns  -------------------------------------

# Set peaking time:  
peMin = 0.02; peMax = 0.033; peakingTime = 140

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); #plt.close()

# Plot PE spectrum:
lower, upper = 0, 140; nbins = 100; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 24; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 160 ns  -------------------------------------

# Set peaking time:  
peMin = 0.02; peMax = 0.03; peakingTime = 180

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 24; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 180 ns  -------------------------------------

# Set peaking time:  
peMin = 0.02; peMax = 0.03; peakingTime = 200

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 15, 24; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 140 ns  -------------------------------------

# Set peaking time:  
peMin = 0.015; peMax = 0.025; peakingTime = 240

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

# Plot PE spectrum:
lower, upper = 0, 100; nbins = 100; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 25; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 200 ns  -------------------------------------

# Set peaking time:  
peMin = 0.014; peMax = 0.023; peakingTime = 280

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 20, 28; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

errors.append(enc_error)
resolutions_2s.append(abs(sig1));
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 220 ns  -------------------------------------

# Set peaking time:  
peMin = 0.012; peMax = 0.021; peakingTime = 320

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 26; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

resolutions_2s.append(abs(sig1));
errors.append(enc_error)
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 270 ns  -------------------------------------

# Set peaking time:  
peMin = 0.01; peMax = 0.018; peakingTime = 400

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 18, 29; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()
#plt.plot(xAxis, gaussianFit, color='red')

delta = []; 
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

resolutions_2s.append(abs(sig1));
errors.append(enc_error)
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 400 ns  -------------------------------------

# Set peaking time:  
peMin = 0.008; peMax = 0.015; peakingTime = 520

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

resolutions_2s.append(abs(sig1));
errors.append(enc_error)
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 500 ns  -------------------------------------

# Set peaking time:  
peMin = 0.006; peMax = 0.013; peakingTime = 680

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

resolutions_2s.append(abs(sig1));
errors.append(enc_error) 
plt.xlim(0, 4);


# --------------------------------- 600 ns  -------------------------------------

# Set peaking time:  
peMin = 0.005; peMax = 0.011; peakingTime = 800

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

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
xMin, xMax = 15, 28; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

resolutions_2s.append(abs(sig1));
errors.append(enc_error)
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 800 ns  -------------------------------------

# Set peaking time:  
peMin = 0.005; peMax = 0.01; peakingTime = 1000

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.04) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 14, 31; 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

resolutions_2s.append(abs(sig1));
errors.append(enc_error)
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 800 ns  -------------------------------------

# Set peaking time:  
peMin = 0.0025; peMax = 0.0064; peakingTime = 1500

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, peMin_unfiltered, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close()

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.04) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 6, 21; 
xR, yR = x[xMin:xMax], y[xMin:xMax]
delt, sig = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0]) 
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
delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax],x,y, expected1); 
mu1 = delt[0]; sig1 = delt[1]; A1 = delt[2]; enc_error = error[1]
print('RESOLUTION: '+str(sig1)+' +/-'+str(enc_error)+'s.p.e.')

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
errors.append(enc_error) 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

plt.close('all')

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs Tpeak plot: 
plt.figure(); 
Tlist_2s = [x*1e-9 for x in Tlist_2s] 
plt.errorbar(Tlist_2s, resolutions_2s, yerr=errors, capthick=0.5, label = 'Data')
#plt.plot(Tlist_2s, resolutions_2s, 'k.', label = 'Data')
plt.plot(tpeak_series_gauss, enc_series_gauss_2sipm, 'k', label = 'Simulations')
plt.title('Gauss filter, 2 SiPMs: 2s, $V_{over}$ = 4.5 V') ; 
plt.xlabel('$t_{peak}$ (s)'); plt.ylabel('enc (s.p.e.)');
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()
plt.savefig('enc-gaussFilter-2s-20k.png')

# Make pandas dataframes: 
df_freq_2s = pd.DataFrame({'peaking time (s)': Tlist_2s,
                   'resolution (s.p.e.)': resolutions_2s,
                   'error (s.p.e.)': errors})
                   
# Save as .csv file: 
df_freq_2s.to_csv('enc-gaussFilter-2s-20k.csv')