"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for PE-spectrum
and plots ENC vs n-sipms for BW and Gauss filter to compare with simulations.  
"""
from sipm_signalProcessing import get_amplitudes, plot_filteredPulses, getAmplitudes_gauss, plot_filteredPulses_gauss
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import pandas as pd
import detect_peaks
import fits
 
# SET PATH: 
# cd Documents/nEXO/sipm_analysis/code

# Read data file from Mietek: 
df_parallel_butt = pd.read_csv('enc_spe_in_nsipms_butt_parallel.csv')
df_series_butt = pd.read_csv('enc_spe_in_nsipms_butt_series.csv')
df_parallel_gauss = pd.read_csv('enc_spe_in_nsipms_gauss_parallel.csv')
df_series_gauss = pd.read_csv('enc_spe_in_nsipms_gauss_series.csv')

# Convert pd data into lists: 
enc_parallel_butt = df_parallel_butt['enc'].tolist()
nsipms_parallel_butt = df_parallel_butt['nsipms'].tolist()
enc_series_butt = df_series_butt['enc'].tolist()
nsipms_series_butt = df_series_butt['nsipms'].tolist()
enc_parallel_gauss = df_parallel_gauss['enc'].tolist()
nsipms_parallel_gauss = df_parallel_gauss['nsipms'].tolist()
enc_series_gauss = df_series_gauss['enc'].tolist()
nsipms_series_gauss = df_series_gauss['nsipms'].tolist()

# Initiate lists for observables:
V_list, nsipms_series, nsipms_parallel = [], [], []
resolutions_bw_series, resolutions_bw_parallel =  [], []
resolutions_gauss_series, resolutions_gauss_parallel = [], []

# Set global variables: 
Npulses = 20000; cutoff = 10e6; bins1 = 200; peak = 'single'; 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----------------------------- Butterworth   ----------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# -------------------------------- SERIES --------------------------------------


# ------------------------------------------------------------------------------
# ------------------------------- 6 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'series';

# Set bias voltage: 
vbias = 68; sipmN = '6'; source = 'led';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append((vbias-59)/2.)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 100; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[2]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; # EDIT 
A1 = y[ind_peak[2]];  # EDIT
# Second PE
mu2 = x[ind_peak[4]]; # EDIT
A2 = y[ind_peak[4]];  # EDIT
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
    delt, sig = fits.singleFitter(x[49:71],y[49:71],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[102:123],y[102:123],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
nsipms_series.append(int(sipmN)); resolutions_bw_series.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.5, 4);
plt.close()

# ------------------------------------------------------------------------------
# -------------------------------- 4 SiPMs -------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 68; sipmN = '4'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append((vbias-59)/2.)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = 0, 120; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
    delt, sig = fits.singleFitter(x[34:51],y[34:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[80:96],y[80:96],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
nsipms_series.append(int(sipmN)); resolutions_bw_series.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
plt.close()

# ------------------------------------------------------------------------------
# -------------------------------- 2 SiPMs -------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 68; sipmN = '2'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append((vbias-59)/2.)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = 0, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[34:51],y[34:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[80:96],y[80:96],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
nsipms_series.append(int(sipmN)); resolutions_bw_series.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
connection = 'series';
plt.close()

# ------------------------------- PARALLEL -------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------- 6 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '6'; source = 'led';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 140; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[2]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; # EDIT 
A1 = y[ind_peak[2]];  # EDIT
# Second PE
mu2 = x[ind_peak[4]]; # EDIT
A2 = y[ind_peak[4]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 6p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[43:75],y[43:75],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[88:118],y[88:118],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
nsipms_parallel.append(int(sipmN)); resolutions_bw_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.5, 4);
plt.close()

# ------------------------------------------------------------------------------
# ------------------------------- 5 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '5'; source = 'led';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[2]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; # EDIT 
A1 = y[ind_peak[2]];  # EDIT
# Second PE
mu2 = x[ind_peak[4]]; # EDIT
A2 = y[ind_peak[4]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 5p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[38:64],y[38:64],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[84:105],y[84:105],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
nsipms_parallel.append(int(sipmN)); resolutions_bw_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.5, 4);
plt.close()

# ------------------------------------------------------------------------------
# ------------------------------- 4 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '4'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[2]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; # EDIT 
A1 = y[ind_peak[2]];  # EDIT
# Second PE
mu2 = x[ind_peak[4]]; # EDIT
A2 = y[ind_peak[4]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 4p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[38:64],y[38:64],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[84:105],y[84:105],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
nsipms_parallel.append(int(sipmN)); resolutions_bw_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
plt.close()

# ------------------------------------------------------------------------------
# ------------------------------- 3 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '3'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
# Plot pulses for debugging: 
#plot_filteredPulses(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[2]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; # EDIT 
A1 = y[ind_peak[2]];  # EDIT
# Second PE
mu2 = x[ind_peak[4]]; # EDIT
A2 = y[ind_peak[4]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[38:56],y[38:56],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[82:100],y[82:100],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
nsipms_parallel.append(int(sipmN)); resolutions_bw_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
plt.close()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------- Gaussian -------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# -------------------------------- SERIES --------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------- 6 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'series';

# Set bias voltage: 
vbias = 68; sipmN = '6'; source = 'led';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append((vbias-59)/2.)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, t_peak = 70); 
# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 100; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[3]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[3]]; # EDIT 
A1 = y[ind_peak[3]];  # EDIT
# Second PE
mu2 = x[ind_peak[5]]; # EDIT
A2 = y[ind_peak[5]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gaussian filter');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[44:65],y[44:65],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[89:109],y[89:109],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gaussian filter, peaking time: 80 ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions_gauss_series.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
plt.close()

# ------------------------------------------------------------------------------
# -------------------------------- 4 SiPMs -------------------------------------
# ------------------------------------------------------------------------------

# Set bias voltage: 
vbias = 68; sipmN = '4'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append((vbias-59)/2.)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, t_peak = 70); 
# Plot pulses for debugging (COMMENT OUT WHEN NOT DEBUGGING): 
#plot_filteredPulses_gauss(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = 0, 120; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gaussian filter');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[31:48],y[31:48],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[73:88],y[73:88],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print(print('Gaussian filter, peaking time: 80 ns'))
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions_gauss_series.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
plt.close()

# ------------------------------------------------------------------------------
# -------------------------------- 2 SiPMs -------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 68; sipmN = '2'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append((vbias-59)/2.)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, t_peak = 70); 
# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = 0, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu2 = 2.; #--------EDIT--------- 
A2 = y[ind_peak[3]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gaussian filter');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[32:45],y[32:45],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[73:85],y[73:85],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gaussian filter, peaking time: 80 ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions_gauss_series.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
plt.close()

# ------------------------------- PARALLEL -------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------- 6 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '6'; source = 'led';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, t_peak = 70); 
# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 140; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[2]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; # EDIT 
A1 = y[ind_peak[2]];  # EDIT
# Second PE
mu2 = x[ind_peak[4]]; # EDIT
A2 = y[ind_peak[4]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Gaussian filter');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 6p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:65],y[33:65],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[78:108],y[78:108],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gaussian filter, peaking time: 80 ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions_gauss_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
plt.close()

# ------------------------------------------------------------------------------
# ------------------------------- 5 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '5'; source = 'led';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, t_peak = 70); 
# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[1]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; # EDIT 
A1 = y[ind_peak[1]];  # EDIT
# Second PE
mu2 = x[ind_peak[3]]; # EDIT
A2 = y[ind_peak[3]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Gaussian filter');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 5p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:50],y[33:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[72:90],y[72:90],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gaussian filter, peaking time: 80 ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions_gauss_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);
plt.close()

# ------------------------------------------------------------------------------
# ------------------------------- 4 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '4'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, t_peak = 70); 
# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[1]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[1]]; # EDIT 
A1 = y[ind_peak[1]];  # EDIT
# Second PE
mu2 = x[ind_peak[3]]; # EDIT
A2 = y[ind_peak[3]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Gaussian filter');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 4p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:57],y[33:57],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[74:91],y[74:91],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gaussian filter, peaking time: 80 ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions_gauss_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
plt.close()

# ------------------------------------------------------------------------------
# ------------------------------- 3 SiPMs  -------------------------------------
# ------------------------------------------------------------------------------

connection = 'parallel';

# Set bias voltage: 
vbias = 34; sipmN = '3'; source = 'dark';
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
V_list.append(vbias-29.5)

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, t_peak = 70); 
# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias) 

# Plot PE spectrum in mV units:
lower, upper = -5, 180; nbins = 200; 
#x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
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
mu1 = x[ind_peak[2]]; # EDIT
x = x/mu1;

# Find estimated values for PE plot fit 
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) 
# Debugging -uncomment plt.close() when done debugging 
plt.figure(); plt.plot(x[0:len(y)],y); 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Single PE
mu1 = x[ind_peak[2]]; # EDIT 
A1 = y[ind_peak[2]];  # EDIT
# Second PE
mu2 = x[ind_peak[4]]; # EDIT
A2 = y[ind_peak[4]];  # EDIT
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str(vbias-29.5)+' V, '+'Gaussian filter');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[35:50],y[35:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[74:90],y[74:90],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gaussian filter, peaking time: 80 ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions_gauss_parallel.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 4);
plt.close()

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs nsipms plot: 
#plt.figure(); 
fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.plot(nsipms_series_butt, enc_series_butt, 'g', label = 'simulations, series' ); #plt.plot(freq_johny, enc_ap, label = 'no afterpulsing' )
plt.plot(nsipms_parallel_butt, enc_parallel_butt, 'k', label = 'simulations, parallel' );
plt.plot(nsipms_series, resolutions_bw_series, 'g.', label = 'data, series' ); #plt.plot(freq_johny, enc_ap, label = 'no afterpulsing' )
plt.plot(nsipms_parallel, resolutions_bw_parallel, 'k.', label = 'data, parallel' );
plt.xlabel('Number of SiPMs'); plt.ylabel('ENC (s.p.e.)'); plt.title('Butterworth filter')
plt.grid(True); plt.show(); plt.legend()

fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.plot(nsipms_series_gauss, enc_series_gauss, 'g', label = 'simulations, series' ); #plt.plot(freq_johny, enc_ap, label = 'no afterpulsing' )
plt.plot(nsipms_parallel_gauss, enc_parallel_gauss, 'k', label = 'simulations, parallel' );
plt.plot(nsipms_series, resolutions_gauss_series, 'g.', label = 'data, series' ); #plt.plot(freq_johny, enc_ap, label = 'no afterpulsing' )
plt.plot(nsipms_parallel, resolutions_gauss_parallel, 'k.', label = 'data, parallel' );
plt.xlabel('Number of SiPMs'); plt.ylabel('ENC (s.p.e.)'); plt.title('Gaussian filter')
plt.grid(True); plt.show(); plt.legend()

"""
# Make pandas dataframes: 
df_freq_noAP = pd.DataFrame({'frequency cut (MHz)': frequencies, 'OV (V)': V_list,
                   'resolution (s.p.e.)': resolutions_noAP})
                   
#df_time = pd.DataFrame({'t_peak (ns)': Tlist,
                   #'resolution (s.p.e.)': resolutions})
                   
df_freq_AP = pd.DataFrame({'frequency cut (MHz)': frequencies, 'OV (V)': V_list,
                   'resolution (s.p.e.)': resolutions})                      
# Save as .csv file: 
df_freq_noAP.to_csv('enc-volt-noAP.csv')
df_freq_AP.to_csv('enc-volt-AP.csv')

#df_time.to_csv('enc-trise-noAP.csv')
"""