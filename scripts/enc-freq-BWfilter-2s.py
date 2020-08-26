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
enc_parallel_butt_2sipm = df_parallel_butt['2sipm_enc'].tolist()

freq_series_butt = df_series_butt['freq'].tolist()
enc_series_butt_2sipm = df_series_butt['2sipm_enc'].tolist()

# Initiate lists for observables:
Tlist_2s, resolutions_2s, frequencies_2s = [], [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; bins1 = 200; peak = 'single'; 


# ------------------------------------------------------------------------------
# ----------------------------- 2-sipms 1s1p  ----------------------------------
# ------------------------------------------------------------------------------

sipmN = '2'; source = 'dark';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# ------------------------------- 800 KHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 0.8e6; peMin = 0.005; peMax = 0.011;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

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

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[0]]; #--------EDIT--------- 
A1 = y[ind_peak[0]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[1]]; #--------EDIT--------- 
A2 = y[ind_peak[1]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[13:26],y[13:26],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[35:47],y[35:47],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 1 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 1e6; peMin = 0.005; peMax = 0.013;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[18:30],y[18:30],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[41:57],y[41:57],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 2 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 2e6; peMin = 0.013; peMax = 0.022;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.01) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.01) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[0]]; #--------EDIT--------- 
A1 = y[ind_peak[0]];  #--------EDIT--------- 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[27:43],y[27:43],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[65:77],y[65:77],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 3 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 3e6; peMin = 0.02; peMax = 0.028;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[40:55],y[40:55],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[89:104],y[89:104],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 4 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 4e6; peMin = 0.02; peMax = 0.033;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.01) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[40:55],y[40:55],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[89:104],y[89:104],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 5 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 5e6; peMin = 0.02; peMax = 0.035;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.02) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:52],y[39:52],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[85:99],y[85:99],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 6 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 6e6; peMin = 0.02; peMax = 0.04;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 160; nbins = 200; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.01) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
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
mu2 = 2; #--------EDIT--------- 
A2 = 0.01;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[37:50],y[37:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[81:92],y[81:92],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 7 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 7e6; peMin = 0.02; peMax = 0.042;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 180; nbins = 200; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.01) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
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
mu2 = 2; #--------EDIT--------- 
A2 = 0.01;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[34:46],y[34:46],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[76:88],y[76:88],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 8 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 8e6; peMin = 0.025; peMax = 0.044;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 180; nbins = 200; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.01) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
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
mu2 = 2; #--------EDIT--------- 
A2 = 0.01;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[35:48],y[35:48],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[80:92],y[80:92],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 9 MHz  --------------------------------------

# Set cutoff frequency: 
cutoff = 9e6; peMin = 0.03; peMax = 0.044;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 180; nbins = 200; 
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.rcParams.update({'font.size': 12})
y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (lower,upper), normed=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
plt.legend(); plt.show()

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15, mph = 0.01) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
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
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 
A2 = 0.01;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[36:49],y[36:49],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[81:94],y[80:94],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 10 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 10e6; peMin = 0.025; peMax = 0.046;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 180; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:50],y[39:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[83:99],y[83:99],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 15 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 15e6; peMin = 0.025; peMax = 0.05;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 180; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:53],y[39:53],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[87:104],y[87:104],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 20 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 20e6; peMin = 0.025; peMax = 0.05;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 180; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:53],y[39:53],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[87:104],y[87:104],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 40 MHz  -------------------------------------

# Set cutoff frequency: 
cutoff = 40e6; peMin = 0.025; peMax = 0.05;

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias, cutoffFreq = cutoff) 

# Plot PE spectrum:
lower, upper = 0, 180; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:52],y[39:52],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[82:98],y[82:98],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies_2s.append(cutoff)   
resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

#plt.figure(); plt.plot(y)
plt.close('all')

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs Tpeak plot: 
plt.figure(); 
plt.plot(frequencies_2s, resolutions_2s, 'k.', label = 'Data')
plt.plot(freq_series_butt, enc_series_butt_2sipm, 'k', label = 'Simulations')
plt.title('BW filter, 2 SiPMs: 2s, $V_{over}$ = 4.5 V'); plt.xlabel('$F_{cut}$ (MHz)'); plt.ylabel('enc (s.p.e.)'); 
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()

# Make ENC vs Tpeak plot: 
plt.figure(); 
plt.plot(Tlist_2s, resolutions_2s, 'g.', label = 'Data')
plt.title('BW filter, 4 SiPMs: 2s2p, $V_{over}$ = 4.5 V'); plt.xlabel('$t_{peak}$ (ns)'); plt.ylabel('enc (s.p.e.)'); 
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()

# Make pandas dataframes: 
df_freq_2s = pd.DataFrame({'frequency cut (MHz)': frequencies_2s, 'peaking time (ns)': Tlist_2s,
                   'resolution (s.p.e.)': resolutions_2s})
                   
# Save as .csv file: 
df_freq_2s.to_csv('enc-bwFilter-measurement-2s.csv')