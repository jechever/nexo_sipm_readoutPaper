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
enc_parallel_butt_6sipm = df_parallel_butt['6sipm_enc'].tolist()

freq_series_butt = df_series_butt['freq'].tolist()
enc_series_butt_6sipm = df_series_butt['6sipm_enc'].tolist()

# Initiate lists for observables:
Tlist_6s, frequencies_6s, resolutions_6s = [], [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; bins1 = 200; peak = 'single'; 

# ------------------------------------------------------------------------------
# ----------------------------- 6-sipms 2s3p  ----------------------------------
# ------------------------------------------------------------------------------
sipmN = '6'; source = 'led';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# ------------------------------- 0.8 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 0.8e6; peMin = 0.003; peMax = 0.007;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 25; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[58:75],y[58:75],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[96:110],y[96:110],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);


# -------------------------------- 1 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 1e6; peMin = 0.003; peMax = 0.008;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 25; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[67:85],y[67:85],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[112:126],y[112:126],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

# -------------------------------- 2 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 2e6; peMin = 0.009; peMax = 0.015;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 50; nbins = 200; 
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
mu1 = 1; #--------EDIT--------- 
A1 = 0.06;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 
A2 = 0.04;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[50:70],y[50:70],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[97:110],y[97:110],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
#plt.figure(); plt.plot(y)

# -------------------------------- 3 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 3e6; peMin = 0.01; peMax = 0.02;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 80; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[41:55],y[41:55],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[78:91],y[78:91],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

# -------------------------------- 4 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 4e6; peMin = 0.01; peMax = 0.024;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 80; nbins = 200; 
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
mu1 = 1; #--------EDIT--------- 
A1 = 0.06;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 
A2 = 0.04;  #--------EDIT---------
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[45:66],y[45:66],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[90:108],y[90:108],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
#plt.figure(); plt.plot(y)

# -------------------------------- 5 MHz  --------------------------------------
# Set cutoff frequency: 
cutoff = 5e6; peMin = 0.01; peMax = 0.025;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 80; nbins = 200; 
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
mu2 = x[ind_peak[5]]; #--------EDIT--------- 
A2 = y[ind_peak[5]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[53:69],y[53:69],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[104:119],y[104:119],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

# -------------------------------- 6 MHz  --------------------------------------
# Set cutoff frequency: 
cutoff = 6e6; peMin = 0.015; peMax = 0.028;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 90; nbins = 200; 
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
mu1 = 1;  
A1 = 0.06;  
# Second PE
mu2 = 2; 
A2 = 0.04;  
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[47:71],y[47:71],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[96:115],y[96:115],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
#plt.figure(); plt.plot(y)

# -------------------------------- 7 MHz  --------------------------------------
# Set cutoff frequency: 
cutoff = 7e6; peMin = 0.015; peMax = 0.03;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, min_rise = 0.02, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 90; nbins = 200; 
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
mu1 = 1;  
A1 = 0.06;  
# Second PE
mu2 = 2; 
A2 = 0.04;  
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[50:75],y[50:75],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[100:122],y[100:122],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
#plt.figure(); plt.plot(y)

# -------------------------------- 8 MHz  --------------------------------------
# Set cutoff frequency: 
cutoff = 8e6; peMin = 0.015; peMax = 0.033;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, min_rise = 0.01, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 100; nbins = 200; 
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
mu1 = 1;  
A1 = 0.06;  
# Second PE
mu2 = 2; 
A2 = 0.04;  
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[45:71],y[45:71],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[94:118],y[94:118],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
#plt.figure(); plt.plot(y)

# -------------------------------- 9 MHz  --------------------------------------
# Set cutoff frequency: 
cutoff = 9e6; peMin = 0.015; peMax = 0.034;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
mu1 = 1;  
A1 = 0.06;  
# Second PE
mu2 = 2; 
A2 = 0.04;  
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:61],y[39:61],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[82:103],y[82:103],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
#plt.figure(); plt.plot(y)

# -------------------------------- 10 MHz  -------------------------------------
# Set cutoff frequency: 
cutoff = 10e6; peMin = 0.015; peMax = 0.035;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:61],y[33:61],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[86:105],y[86:105],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

# -------------------------------- 15 MHz  -------------------------------------
# Set cutoff frequency: 
cutoff = 15e6; peMin = 0.015; peMax = 0.037;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, min_rise = 0.03, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[45:65],y[45:65],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[88:109],y[88:109],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 20 MHz  -------------------------------------
# Set cutoff frequency: 
cutoff = 20e6; peMin = 0.02; peMax = 0.038;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, min_rise = 0.03, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[43:65],y[43:65],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[88:110],y[88:110],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 40 MHz  -------------------------------------
# Set cutoff frequency: 
cutoff = 40e6; peMin = 0.02; peMax = 0.038;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, min_rise = 0.03, cutoffFreq = cutoff, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[41:65],y[41:65],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[81:106],y[81:106],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_6s.append(cutoff)   
resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

plt.close('all')

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs Tpeak plot: 
plt.figure(); 
#plt.plot(frequencies, resolutions, 'k.')
plt.plot(frequencies_6s, resolutions_6s, 'b.', label = 'Data')
plt.plot(freq_series_butt, enc_series_butt_6sipm, 'b', label = 'Simulations')
plt.title('BW filter, 6 SiPMs: 2s3p, $V_{over}$ = 4.5 V') ; plt.xlabel('$F_{cut}$ (MHz)'); plt.ylabel('enc (s.p.e.)');
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()

# Make ENC vs Tpeak plot: 
plt.figure(); 
plt.plot(Tlist_6s, resolutions_6s, 'g.', label = 'Data')
plt.title('BW filter, 4 SiPMs: 2s2p, $V_{over}$ = 4.5 V'); plt.xlabel('$t_{peak}$ (ns)'); plt.ylabel('enc (s.p.e.)'); 
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()

# Make pandas dataframes: 
df_freq_6s = pd.DataFrame({'frequency cut (MHz)': frequencies_6s, 'peaking time (ns)': Tlist_6s,
                   'resolution (s.p.e.)': resolutions_6s})
                   
# Save as .csv file: 
df_freq_6s.to_csv('enc-bwFilter-measurement-6s.csv')