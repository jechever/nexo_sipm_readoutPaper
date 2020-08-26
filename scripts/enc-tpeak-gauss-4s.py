"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 4-SiPMs 
(2s2p) for t-rise dependence study, using Gaussian filter. Uses sample biased 
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
#df_parallel_gauss = pd.read_csv('simulations_mietek/enc_semigauss_shaper_parallel.csv')
df_series_gauss = pd.read_csv('simulations_mietek/enc_semigauss_shaper_2_in_series.csv')

tpeak_series_gauss = df_series_gauss['tpeak'].tolist()
enc_series_gauss_4sipm = df_series_gauss['4sipm_enc'].tolist()

# Initiate lists for observables:
resolutions_4s, Tlist_4s = [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; bins1 = 200; peak = 'single';

# ------------------------------------------------------------------------------
# ----------------------------- 4-sipms 2s2p  ----------------------------------
# ------------------------------------------------------------------------------
sipmN = '4'; source = 'dark';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# --------------------------------- 50 ns  -------------------------------------

# Set peaking time:  
peMin = 0.017; peMax = 0.032; peakingTime = 50

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 
plt.close()

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[40:57],y[40:57],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[79:98],y[79:98],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));                  
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

#plt.figure(); plt.plot(y)

# --------------------------------- 80 ns  -------------------------------------

# Set peaking time:  
peMin = 0.017; peMax = 0.028; peakingTime = 80

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[35:51],y[35:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[75:91],y[75:91],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[33:50],y[33:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[69:84],y[69:84],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[30:47],y[30:47],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[67:80],y[67:80],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[29:45],y[29:45],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[62:77],y[62:77],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[29:42],y[29:42],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[61:72],y[61:72],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[28:41],y[29:41],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[56:68],y[56:68],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 120; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[26:39],y[26:39],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[54:64],y[54:64],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 80; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[36:53],y[36:53],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[74:87],y[74:87],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 60; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[45:64],y[45:64],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[89:104],y[89:104],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 60; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[43:61],y[43:61],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[82:97],y[82:97],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
    delt, sig = fits.singleFitter(x[36:63],y[36:63],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[96:112],y[96:112],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# --------------------------------- 400 ns  -------------------------------------

# Set peaking time:  
peMin = 0.005; peMax = 0.01; peakingTime = 520

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
    delt, sig = fits.singleFitter(x[27:53],y[27:53],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[70:97],y[70:97],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
    delt, sig = fits.singleFitter(x[20:44],y[20:44],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[53:78],y[53:78],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[21:52],y[21:52],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[64:89],y[64:89],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 200; 
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
mu1 = 1; 
A1 = 0.1;   
# Second PE
mu2 = 2; 
A2 = 0.05;  
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
    delt, sig = fits.singleFitter(x[17:47],y[17:47],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[56:79],y[56:79],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 330 ns  -------------------------------------

# Set peaking time:  
peMin = 0.001; peMax = 0.006; peakingTime = 1000

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 200; 
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
mu1 = 1; 
A1 = 0.1;   
# Second PE
mu2 = 2; 
A2 = 0.05;  
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
    delt, sig = fits.singleFitter(x[16:43],y[16:43],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[50:71],y[50:71],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
#plt.figure(); plt.plot(y)

# --------------------------------- 1000 ns  -------------------------------------

# Set peaking time:  
peMin = 0.001; peMax = 0.005; peakingTime = 1300

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.01, Header = header); Tlist_4s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = 0, 30; nbins = 200; 
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
mu1 = 1; 
A1 = 0.1;   
# Second PE
mu2 = 2; 
A2 = 0.05;  
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
    delt, sig = fits.singleFitter(x[9:36],y[9:36],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[41:56],y[41:56],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_4s.append(abs(sig1));
#resolutions_4s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);
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

# Make pandas dataframes: 
df_freq_4s = pd.DataFrame({'peaking time (s)': Tlist_4s,
                   'resolution (s.p.e.)': resolutions_4s})
                   
# Save as .csv file: 
df_freq_4s.to_csv('enc-gaussFilter-measurement-4s.csv')