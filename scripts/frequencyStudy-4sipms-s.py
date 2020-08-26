"""
Full analysis script for 2019 SiPM paper 
Makes resolution fits for 4-SiPM series connection for freq. dependence study
"""
from sipm_signalProcessing import get_amplitudes, plot_peSpectrum, get_amplitudesFullWindow
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import sipm_dataAnalysis
import detect_peaks
import pylab
import fits
 
# SET PATH: 
# cd Documents/nEXO/sipm_analysis/code
# Dark data only 


# ------------------------------------------------------------------------------
# --------------------------------- 1 MHz  -------------------------------------
# ------------------------------------------------------------------------------

Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 1e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum:
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys #2 
hist_lowerBound = ind_valleys[0]; #EDIT
hist_upperBound = ind_valleys[5]  #EDIT
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
#pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

"""
if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[10:19],y[10:19],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[33:42],y[33:42],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2))) 
"""   
plt.xlim(0, 5); plt.ylim(-5, 1800) 

frequencies = []; resolutions = []; 
# ------------------------------------------------------------------------------
# --------------------------------- 3 MHz  -------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 3e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 

# Plot PE spectrum:
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #EDIT
hist_upperBound = ind_valleys[5]  #EDIT
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[10:19],y[10:19],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[33:42],y[33:42],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    error = np.sqrt(); 
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2))+'+/- '+str(error))
    
frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2)); 
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# --------------------------------- 5 MHz  -------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 5e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #--------EDIT--------- 
hist_upperBound = ind_valleys[6]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[7:16],y[7:16],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[30:38],y[30:38],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str())
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));    
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# -------------------------------- 10 MHz  -------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 10e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #--------EDIT--------- 
hist_upperBound = ind_valleys[9]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[4:18],y[4:18],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[30:41],y[30:41],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));    
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# -------------------------------- 20 MHz  -------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 20e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
# ---------------------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #--------EDIT--------- 
hist_upperBound = ind_valleys[7]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$= -'+str(vbias)+' V, '+'Filter= '+str(int(cutoff/(10**6)))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[10:28],y[10:28],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[42:54],y[42:54],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));
plt.xlim(0, 5); plt.ylim(-5, 1800)  

# ------------------------------------------------------------------------------
# -------------------------------- 40 MHz  -------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 40e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
# ---------------------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #--------EDIT--------- 
hist_upperBound = ind_valleys[7]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$= -'+str(vbias)+' V, '+'Filter= '+str(int(cutoff/(10**6)))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[3:13],y[3:13],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[26:39],y[26:39],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# -------------------------------- 60 MHz  -------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 60e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #--------EDIT--------- 
hist_upperBound = ind_valleys[9]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[12:24],y[12:24],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[35:49],y[35:49],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# -------------------------------- 80 MHz  -------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 80e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #--------EDIT--------- 
hist_upperBound = ind_valleys[6]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS

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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[18:30],y[18:30],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[45:58],y[45:58],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))
    
frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));    
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# -------------------------------- 100 MHz  ------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 100e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #--------EDIT--------- 
hist_upperBound = ind_valleys[6]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
# ------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS

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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[7:18],y[7:18],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[33:46],y[33:46],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2)))

frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));     
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# -------------------------------- 120 MHz  ------------------------------------
# ------------------------------------------------------------------------------
Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
bins1 = 120; peak = 'single'; cutoff = 120e6; 
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, cutoffFreq = cutoff); 
#plot_peSpectrum
lower, upper = 0, 200; nbins = 200; 
x, y = plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);

# Find peaks and valleys #1 
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[0]; #EDIT
hist_upperBound = ind_valleys[9]  #EDIT
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE4 = x[ind_peak[5]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(amplitudes, bins1, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
plt.close()

# EDIT X-AXIS
# Normalize by dividing by first PE peak value to convert mV -> PE 
mu1 = x[ind_peak[0]];
x = x/mu1;
# ------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
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
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs in '+connection+', V$_{bias}$ = -'+str(vbias)+' V, '+'Filter: '+str(cutoff/(10**6))+' MHz');
plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.legend(loc='upper right');  plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[14:23],y[14:23],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[35:47],y[35:47],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('Filter: '+str(cutoff/(10**6))+' MHz')
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2))) 

frequencies.append(cutoff/(10**6)); resolutions.append(abs(sig1)/abs(mu1-mu2));
plt.xlim(0, 5); plt.ylim(-5, 1800) 

# ------------------------------------------------------------------------------
# ---------------------------- Resolution plot ---------------------------------
# ------------------------------------------------------------------------------

import pandas as pd 

# Make ENC vs Freq. plot: 
plt.figure(); plt.plot(frequencies, resolutions, 'k.')
plt.xlabel('Cutoff frequency (MHz)'); plt.ylabel('ENC (s.p.e.)'); 
plt.grid(); plt.show()

# Make pandas dataframe: 
df = pd.DataFrame({'frequency cut (MHz)': frequencies,
                   'resolution (s.p.e.)': resolutions})
# Save as .csv file: 
df.to_csv('enc-freq.csv', sep='\t')

# Save file: 