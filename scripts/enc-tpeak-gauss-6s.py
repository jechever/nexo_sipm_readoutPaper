"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 6-SiPMs 
(2s3p) for t-rise dependence study, using Gaussian filter. Uses sample biased 
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
enc_series_gauss_6sipm = df_series_gauss['6sipm_enc'].tolist()

# Initiate lists for observables:
resolutions_6s, Tlist_6s  = [], []

# Set global variables: 
Npulses = 20000; vbias = 68; connection = 'series'; bins1 = 200; peak = 'single';

# ------------------------------------------------------------------------------
# ----------------------------- 6-sipms 2s3p  ----------------------------------
# ------------------------------------------------------------------------------
sipmN = '6'; source = 'led';
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

# --------------------------------- 50 ns  -------------------------------------

# Set peaking time:  
peMin = 0.015; peMax = 0.034; peakingTime = 50

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.02, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[40:72],y[40:72],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[87:118],y[87:118],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# -------------------------------- 60 ns  -------------------------------------

# Set peaking time:  
peMin = 0.015; peMax = 0.028; peakingTime = 80

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, min_rise = 0.02, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[40:64],y[40:64],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[83:108],y[83:108],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[40:59],y[40:59],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[83:100],y[83:100],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[37:55],y[37:55],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[78:92],y[78:92],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[37:54],y[37:54],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[74:86],y[74:86],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[35:51],y[35:51],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[70:82],y[70:82],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[33:48],y[33:48],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[66:78],y[66:78],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[30:47],y[30:47],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[63:74],y[63:74],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[29:43],y[29:43],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[58:68],y[58:68],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[2]];
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[44:63],y[44:63],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[85:100],y[85:100],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
ind_peak = detect_peaks.detect_peaks(y, mpd = 15) # peaks in PHS
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[42:59],y[42:59],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[78:94],y[78:94],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 40; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[55:76],y[55:76],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[100:118],y[100:118],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 40; nbins = 200; 
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[47:66],y[47:66],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[86:101],y[86:101],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 30; nbins = 200; 
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
mu1 = x[ind_peak[3]];
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[54:74],y[54:74],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[90:109],y[90:109],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
mu1 = x[ind_peak[3]];
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[59:80],y[59:80],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[98:117],y[98:117],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[56:76],y[56:76],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[90:106],y[90:106],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

#plt.figure(); plt.plot(y)

# ------------------------------- 1000 ns  -------------------------------------

# Set peaking time:  
peMin = 0.0025; peMax = 0.0065; peakingTime = 1000

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
pulse, Trise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, peakingTime, Header = header); Tlist_6s.append(Trise)
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[52:70],y[52:70],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[83:98],y[83:98],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 16; nbins = 200; 
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
mu1 = x[ind_peak[3]];
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
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 3p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[71:93],y[71:93],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[106:128],y[106:128],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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
amplitudes, integrals = getAmplitudes_gauss(pathToData, Npulses, header, vbias, peakingTime, fullWindow = True); 

# Plot pulses for debugging: 
#plot_filteredPulses_gauss(pathToData, header, vbias, t_peak = Trise) 

# Plot PE spectrum:
lower, upper = -5, 20; nbins = 200; 
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
mu1 = x[ind_peak[3]];
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
expected2 = [mu2,sig,A2];# Expected fit values 

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
    delt, sig = fits.singleFitter(x[55:70],y[55:70],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[80:92],y[80:92],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Gauss filter: '+str(Trise)+' ns')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions_6s.append(abs(sig1));
#resolutions_6s.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(-0.2, 4);

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

# Make pandas dataframes: 
df_freq_6s = pd.DataFrame({'peaking time (ns)': Tlist_6s,
                   'resolution (s.p.e.)': resolutions_6s})
                   
# Save as .csv file: 
df_freq_6s.to_csv('enc-gaussFilter-measurement-6s.csv')