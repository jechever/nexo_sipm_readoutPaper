"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 6-SiPMs 
(2s3p) for t-rise dependence study, using Gaussian filter. Uses sample biased 
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
enc_series_gauss_6sipm = df_series_gauss['6sipm_enc'].tolist()

# Initiate lists for observables:
resolutions_6s, Tlist_6s  = [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; bins1 = 200; peak = 'single';
peMin_unfiltered = 0.02; 

# ------------------------------------------------------------------------------
# ----------------------------- 6-sipms 2s3p  ----------------------------------
# ------------------------------------------------------------------------------
sipmN = '6'; source = 'led';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
"""
# --------------------------------- 50 ns  -------------------------------------

# Set peaking time:  
peMin = 0.015; peMax = 0.034; peakingTime = 50

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.02, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[16:31],y[16:31],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    #delt, sig = fits.singleFitter(x[87:118],y[87:118],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    print('RESOLUTION: '+str(sig1))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)
"""
# -------------------------------- 60 ns  -------------------------------------

# Set peaking time:  
peMin = 0.015; peMax = 0.028; peakingTime = 80

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.02, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise); plt.close() 

# Plot PE spectrum:
lower, upper = -5, 100; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 18, 33; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 80 ns  -------------------------------------

# Set peaking time:  
peMin = 0.015; peMax = 0.028; peakingTime = 100

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 100; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 17, 31; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 100 ns  -------------------------------------

# Set peaking time:  
peMin = 0.012; peMax = 0.026; peakingTime = 120

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 100; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 17, 29; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 120 ns  -------------------------------------

# Set peaking time:  
peMin = 0.012; peMax = 0.025; peakingTime = 140

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 100; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 17, 26; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 140 ns  -------------------------------------

# Set peaking time:  
peMin = 0.012; peMax = 0.025; peakingTime = 160

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 100; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 16, 25; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 160 ns  -------------------------------------

# Set peaking time:  
peMin = 0.009; peMax = 0.022; peakingTime = 180

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 100; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 15, 24; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 170 ns  -------------------------------------

# Set peaking time:  
peMin = 0.008; peMax = 0.020; peakingTime = 200

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 13, 24; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 200 ns  -------------------------------------

# Set peaking time:  
peMin = 0.008; peMax = 0.018; peakingTime = 240

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 12, 21; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 220 ns  -------------------------------------

# Set peaking time:  
peMin = 0.007; peMax = 0.017; peakingTime = 280

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 14, 25; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 260 ns  -------------------------------------

# Set peaking time:  
peMin = 0.007; peMax = 0.015; peakingTime = 320

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 13, 23; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 300 ns  -------------------------------------

# Set peaking time:  
peMin = 0.006; peMax = 0.013; peakingTime = 400

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 17, 28; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 400 ns  -------------------------------------

# Set peaking time:  
peMin = 0.005; peMax = 0.010; peakingTime = 520

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 13, 23; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 500 ns  -------------------------------------

# Set peaking time:  
peMin = 0.004; peMax = 0.008; peakingTime = 680

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 13, 26; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 600 ns  -------------------------------------

# Set peaking time:  
peMin = 0.0035; peMax = 0.007; peakingTime = 780

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 13, 28; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 720 ns  -------------------------------------

# Set peaking time:  
peMin = 0.0035; peMax = 0.0065; peakingTime = 900

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 10, 25; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# ------------------------------- 820 ns  -------------------------------------

# Set peaking time:  
peMin = 0.0025; peMax = 0.0065; peakingTime = 1100

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 8, 21; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# ------------------------------- 1000 ns  -------------------------------------

# Set peaking time:  
peMin = 0.0023; peMax = 0.005; peakingTime = 1250

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes_noAP(pathToData, Npulses, header, vbias, peMin_unfiltered, t_peak = peakingTime, fullWindow = True);
print('Number of waveforms: '+str(len(amplitudes)))

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 16; nbins = 100; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 10, 31; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)

# ------------------------------- 1500 ns  -------------------------------------

# Set peaking time:  
peMin = 0.0017; peMax = 0.0046; peakingTime = 1500

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

# Make initial fit: 
delta = []; 
# RANGES: 
xMin, xMax = 4, 20; 
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

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);

#plt.figure(); plt.plot(y)
plt.close('all')

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs Tpeak plot: 
plt.figure(); 
Tlist_6s = [x*1e-9 for x in Tlist_6s] 
plt.plot(Tlist_6s, resolutions_6s, 'b.', label = 'Data')
plt.plot(tpeak_series_gauss, enc_series_gauss_6sipm, 'b', label = 'Simulations')
plt.title('Gauss filter, 6 SiPMs: 2s3p, $V_{over}$ = 4.5 V') ; 
plt.xlabel('$t_{peak}$ (s)'); plt.ylabel('enc (s.p.e.)');
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()
plt.savefig('enc-gaussFilter-noAP-6s.png')

# Make pandas dataframes: 
df_freq_6s = pd.DataFrame({'peaking time (ns)': Tlist_6s,
                   'resolution (s.p.e.)': resolutions_6s})
                   
# Save as .csv file: 
df_freq_6s.to_csv('enc-gaussFilter-noAP-6s.csv')