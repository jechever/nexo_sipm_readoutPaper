"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 2-SiPMs 
(2s) for t-rise dependence study, using Gaussian filter. Uses sample biased 
at Vbias = -68 V (Vover = 4.5 V)
"""
from sipm_signalProcessing import getAvgPulse_gauss, getAmplitudes_gauss, plot_filteredPulses_gauss
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import pandas as pd 
import detect_peaks
import fits
 
# SET PATH: 
# cd Documents/nEXO/sipm_analysis/code

# Read data file from Mietek: 
df_series_gauss = pd.read_csv('simulations_mietek/enc_semigauss_shaper_2_in_series.csv')
tpeak_series_gauss = df_series_gauss['tpeak'].tolist()
enc_series_gauss_2sipm = df_series_gauss['2sipm_enc'].tolist()

# Initiate lists for observables:
resolutions_2s, Tlist_2s = [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; bins1 = 200; peak = 'single';

# ------------------------------------------------------------------------------
# ----------------------------- 2-sipms 1s1p  ----------------------------------
# ------------------------------------------------------------------------------
sipmN = '2'; source = 'dark';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# --------------------------------- 50 ns  -------------------------------------

# Set peaking time:  
peMin = 0.028; peMax = 0.043; peakingTime = 50

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 
plt.close()

# Plot PE spectrum:
lower, upper = 0, 160; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[38:54],y[38:54],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[86:101],y[86:101],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1)); 
plt.xlim(0.2, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 80 ns  -------------------------------------

# Set peaking time:  
peMin = 0.026; peMax = 0.04; peakingTime = 80

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 
plt.close()

# Plot PE spectrum:
lower, upper = 0, 160; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[35:49],y[35:49],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[79:92],y[79:92],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 
plt.close()

# Plot PE spectrum:
lower, upper = 0, 160; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[34:46],y[34:46],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[75:87],y[75:87],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 
plt.close()

# Plot PE spectrum:
lower, upper = 0, 160; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[30:43],y[30:43],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[70:82],y[70:82],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 
plt.close()

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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:47],y[33:47],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[75:89],y[75:89],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 
plt.close()

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[34:49],y[34:49],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[77:92],y[77:92],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 120; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[32:47],y[32:47],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[73:87],y[73:87],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:50],y[33:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[78:94],y[78:94],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 80; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[39:57],y[39:57],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[88:103],y[88:103],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1));
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 80; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:53],y[33:53],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[79:96],y[79:96],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_2s.append(abs(sig1));
#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 60; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[37:60],y[37:60],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[89:111],y[89:111],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_2s.append(abs(sig1));
#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 60; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[28:51],y[28:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[69:89],y[69:89],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_2s.append(abs(sig1));
#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 60; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[21:40],y[21:40],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[57:70],y[57:70],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_2s.append(abs(sig1));
#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);


# --------------------------------- 600 ns  -------------------------------------

# Set peaking time:  
peMin = 0.005; peMax = 0.011; peakingTime = 800

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_2s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 40; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[26:58],y[26:58],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[71:95],y[71:95],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_2s.append(abs(sig1));
#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 200; 
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

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[28:65],y[28:65],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[78:105],y[78:105],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_2s.append(abs(sig1));
#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) ; plt.close()

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 200; 
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
mu1 = x[ind_peak[1]];
x = x/mu1;

# Single PE
mu1 = 1; #--------EDIT--------- 
A1 = 0.05;  #--------EDIT--------- 
# Second PE
mu2 = 2; #--------EDIT--------- 0.0
A2 = 0.02;  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(Trise));
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[9:44],y[9:44],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    #delt, sig = fits.singleFitter(x[64:80],y[64:80],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    print('RESOLUTION: '+str(sig1))

#resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
resolutions_2s.append(abs(sig1)); 
plt.xlim(0.2, 4);

#plt.figure(); plt.plot(y)

plt.close('all')

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

# Make ENC vs Tpeak plot: 
plt.figure(); 
Tlist_4s = [x*1e-9 for x in Tlist_2s] 
plt.plot(Tlist_4s, resolutions_2s, 'k.', label = 'Data')
plt.plot(tpeak_series_gauss, enc_series_gauss_2sipm, 'k', label = 'Simulations')
plt.title('Gauss filter, 2 SiPMs: 2s, $V_{over}$ = 4.5 V') ; 
plt.xlabel('$t_{peak}$ (s)'); plt.ylabel('enc (s.p.e.)');
plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()

# Make pandas dataframes: 
df_freq_2s = pd.DataFrame({'peaking time (s)': Tlist_2s,
                   'resolution (s.p.e.)': resolutions_2s})
                   
# Save as .csv file: 
df_freq_2s.to_csv('enc-gaussFilter-measurement-2s.csv')