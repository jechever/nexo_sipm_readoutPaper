"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 4-SiPMs 
in series connection for freq. and t-rise dependence study. Uses sample biased 
at Vbias = -70 V 
"""
from sipm_signalProcessing import getAvgPulse, get_amplitudes, plot_filteredPulses
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import pandas as pd 
import detect_peaks
import fits
 
# SET PATH: 
# cd Documents/nEXO/sipm_analysis/code


# Read data file from Mietek: 
df_parallel_butt = pd.read_csv('simulations_mietek/enc_spe_in_freq_butt_parallel.csv')
df_series_butt = pd.read_csv('simulations_mietek/enc_spe_in_freq_butt_2_in_series.csv')

# Convert pd data into lists: 
freq_parallel_butt = df_parallel_butt['freq'].tolist()
enc_parallel_butt_4sipm = df_parallel_butt['4sipm_enc'].tolist()

freq_series_butt = df_series_butt['freq'].tolist()
enc_series_butt_4sipm = df_series_butt['4sipm_enc'].tolist()

# Initiate lists for observables:
Tlist_4s, frequencies_4s, resolutions_4s = [], [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; bins1 = 200; peak = 'single'; 


# ------------------------------------------------------------------------------
# ----------------------------- 4-sipms 2s2p  ----------------------------------
# ------------------------------------------------------------------------------

sipmN = '4'; source = 'dark';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# ------------------------------- 800 KHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 0.8e6; peMin = 0.002; peMax = 0.085;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 40; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[2]]; #--------EDIT--------- 
A2 = y[ind_peak[2]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[16:36],y[16:36],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[46:63],y[46:63],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)   
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 1 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 1e6; peMin = 0.003; peMax = 0.01;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 40; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[21:43],y[21:43],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[60:76],y[60:76],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies_4s.append(cutoff)   
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 2 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 2e6; peMin = 0.01; peMax = 0.016;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 60; nbins = 200; 
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
mu1 = x[ind_peak[2]];
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; #--------EDIT--------- 
A1 = y[ind_peak[2]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[4]]; #--------EDIT--------- 
A2 = y[ind_peak[4]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[30:50],y[30:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[76:87],y[76:87],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)   
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 3 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 3e6; peMin = 0.013; peMax = 0.02;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 100; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[24:39],y[24:39],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[60:69],y[60:69],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
  
frequencies_4s.append(cutoff)   
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 4 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 4e6; peMin = 0.015; peMax = 0.023;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 100; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[29:47],y[29:47],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[72:83],y[72:83],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)   
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 5 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 5e6; peMin = 0.018; peMax = 0.026;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 100; nbins = 200; 
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
mu1 = x[ind_peak[2]];
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; #--------EDIT--------- 
A1 = y[ind_peak[2]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[4]]; #--------EDIT--------- 
A2 = y[ind_peak[4]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[35:51],y[35:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[81:94],y[81:94],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies_4s.append(cutoff)     
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 6 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 6e6; peMin = 0.019; peMax = 0.028;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 100; nbins = 200; 
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
mu1 = x[ind_peak[2]];
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; #--------EDIT--------- 
A1 = y[ind_peak[2]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[4]]; #--------EDIT--------- 
A2 = y[ind_peak[4]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[36:56],y[36:56],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[85:102],y[85:102],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)     
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 7 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 7e6; peMin = 0.02; peMax = 0.03;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 100; nbins = 200; 
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
mu1 = x[ind_peak[2]];
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; #--------EDIT--------- 
A1 = y[ind_peak[2]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[4]]; #--------EDIT--------- 
A2 = y[ind_peak[4]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:59],y[39:59],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[92:107],y[92:107],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)     
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 8 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 8e6; peMin = 0.021; peMax = 0.031;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:50],y[33:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[78:94],y[78:94],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)     
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 9 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 9e6; peMin = 0.022; peMax = 0.032;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:54],y[33:54],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[80:99],y[80:99],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)     
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 10 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 10e6; peMin = 0.022; peMax = 0.034;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 140; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[30:48],y[30:48],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[72:88],y[72:88],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
  
frequencies_4s.append(cutoff)   
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 15 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 15e6; peMin = 0.023; peMax = 0.038;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 140; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[30:51],y[30:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[72:93],y[72:93],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)    
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 20 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 20e6; peMin = 0.022; peMax = 0.038;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 140; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[30:51],y[30:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[72:93],y[72:93],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)     
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 40 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 40e6; peMin = 0.022; peMax = 0.037;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, min_rise = 0.03, cutoffFreq = cutoff, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 140; nbins = 200; 
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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; #--------EDIT--------- 
A1 = y[ind_peak[1]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[3]]; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[30:50],y[30:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[72:90],y[72:90],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_4s.append(cutoff)    
resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
plt.close('all')

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs Freq plot: 
plt.figure(); 
plt.plot(frequencies_4s, resolutions_4s, 'g.', label = 'Data')
plt.plot(freq_series_butt, enc_series_butt_4sipm, 'g', label = 'Simulations')
plt.title('BW filter, 4 SiPMs: 2s2p, $V_{over}$ = 4.5 V'); plt.xlabel('$F_{cut}$ (MHz)'); plt.ylabel('enc (s.p.e.)'); 
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()

# Make ENC vs Tpeak plot: 
plt.figure(); 
plt.plot(Tlist_4s, resolutions_4s, 'g.', label = 'Data')
plt.title('BW filter, 4 SiPMs: 2s2p, $V_{over}$ = 4.5 V'); plt.xlabel('$t_{peak}$ (ns)'); plt.ylabel('enc (s.p.e.)'); 
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()

# Make pandas dataframes: 
df_freq_4s = pd.DataFrame({'frequency cut (MHz)': frequencies_4s, 'peaking time (ns)': Tlist_4s,
                   'resolution (s.p.e.)': resolutions_4s})
             
# Save as .csv file: 
df_freq_4s.to_csv('enc-bwFilter-measurement-4s.csv')