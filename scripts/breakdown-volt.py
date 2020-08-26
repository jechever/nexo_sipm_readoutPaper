"""
Full analysis script for 2019 SiPM paper. Makes resolution fits for 4-SiPMs 
in series connection for freq. and t-rise dependence study. Uses sample biased 
at Vbias = -70 V 
"""
from sipm_signalProcessing import getAvgPulse, getAvgPulse2, get_amplitudes_FW_noAP, get_amplitudes_noAP
from sipm_signalProcessing import get_amplitudes, get_amplitudesFullWindow
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import detect_peaks
import fits
 
# SET PATH: 
# cd Documents/nEXO/sipm_analysis/code

# Initiate lists for observables:
V_list, frequencies, resolutions_noAP, resolutions = [], [], [], []

# Set global variables: 
Npulses = 20000; sipmN = '4'; cutoff = 5e6;
source = 'dark'; connection = 'series'; 
bins1 = 200; peak = 'single'; 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------------------- No Afterpulsing  --------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --------------------------------- 70 V  --------------------------------------
# ------------------------------------------------------------------------------

# Set bias voltage: 
vbias = 70; peMin = 0.008; peMax = 0.06;
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
amplitudes, integrals = get_amplitudes_FW_noAP(pathToData, Npulses, header, vbias, peMin, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[43:61],y[43:61],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[98:114],y[98:114],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions_noAP.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);


# ------------------------------------------------------------------------------
# --------------------------------- 69 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 69;  peMin = 0.008; peMax = 0.06;
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
amplitudes, integrals = get_amplitudes_FW_noAP(pathToData, Npulses, header, vbias, peMin, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[40:55],y[40:55],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[90:99],y[90:99],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions_noAP.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);


# ------------------------------------------------------------------------------
# --------------------------------- 68 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 68;  peMin = 0.008; peMax = 0.06;
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
amplitudes, integrals = get_amplitudes_FW_noAP(pathToData, Npulses, header, vbias, peMin, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[34:50],y[34:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[80:95],y[80:95],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions_noAP.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# ------------------------------------------------------------------------------
# --------------------------------- 67 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 67;  peMin = 0.008; peMax = 0.06;
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
amplitudes, integrals = get_amplitudes_FW_noAP(pathToData, Npulses, header, vbias, peMin, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[25:45],y[25:45],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[70:85],y[70:85],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions_noAP.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# ------------------------------------------------------------------------------
# --------------------------------- 66 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 66;  peMin = 0.008; peMax = 0.06;
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
amplitudes, integrals = get_amplitudes_FW_noAP(pathToData, Npulses, header, vbias, peMin, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[25:41],y[25:41],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[62:75],y[62:75],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions_noAP.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# ------------------------------------------------------------------------------
# --------------------------------- 65 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 65;  peMin = 0.008; peMax = 0.06;
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
amplitudes, integrals = get_amplitudes_FW_noAP(pathToData, Npulses, header, vbias, peMin, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[21:38],y[21:38],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[54:66],y[54:66],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions_noAP.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------------------------- With Afterpulsing  -------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# --------------------------------- 70 V  --------------------------------------
# ------------------------------------------------------------------------------

# Set bias voltage: 
vbias = 70; peMin = 0.008; peMax = 0.06;
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[43:61],y[43:61],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[98:114],y[98:114],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

resolutions.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);


# ------------------------------------------------------------------------------
# --------------------------------- 69 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 69;  peMin = 0.008; peMax = 0.06;
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[38:57],y[38:57],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[90:105],y[90:105],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);


# ------------------------------------------------------------------------------
# --------------------------------- 68 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 68;  peMin = 0.008; peMax = 0.06;
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[34:50],y[34:50],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[80:95],y[80:95],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# ------------------------------------------------------------------------------
# --------------------------------- 67 V  --------------------------------------
# ------------------------------------------------------------------------------

# Set bias voltage: 
vbias = 67;  peMin = 0.008; peMax = 0.06;
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz BW');
plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2p2s')
plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[25:45],y[25:45],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[70:85],y[70:85],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# ------------------------------------------------------------------------------
# --------------------------------- 66 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 66;  peMin = 0.008; peMax = 0.06;
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[25:41],y[25:41],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[62:75],y[62:75],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);


# ------------------------------------------------------------------------------
# --------------------------------- 65 V  --------------------------------------
# ------------------------------------------------------------------------------


# Set bias voltage: 
vbias = 65;  peMin = 0.008; peMax = 0.06;
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Path to dataset: 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

if source == 'dark': 
    header = 1;
else: 
    header = 0; 

# Compute rise time and get amplitudes 
#pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias, cutoffFreq = cutoff); Tlist.append(Trise)
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum in mV units:
lower, upper = 0, 100; nbins = 200; 
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
    delt, sig = fits.singleFitter(x[21:38],y[21:38],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[54:66],y[54:66],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    #error = np.sqrt(); 
    #print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
resolutions.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0.5, 4);

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

import pandas as pd 

# Make ENC vs Tpeak plot: 
plt.figure(); 
#plt.plot(frequencies, resolutions, 'k.')
plt.plot(V_list, resolutions, 'r.', label = 'Data with afterpulsing')
plt.plot(V_list, resolutions_noAP, 'k.', label = 'Data with no afterpulsing')
plt.xlabel('$V_{over}$ (V)'); plt.ylabel('enc (s.p.e.)'); 
plt.grid(); plt.show(); plt.legend()

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
