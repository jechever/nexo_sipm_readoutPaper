"""
Full analysis script for 2019 SiPM paper 
Makes resolution fits for 4-SiPM series connection
"""
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import sipm_dataAnalysis
import detect_peaks
import pylab
import fits
 
# ------------------------------------------
# -------- ALWAYS SET PATH FIRST  ----------
# ------------------------------------------

#  cd Documents/nEXO/sipm_analysis/code


"""
# ------------------------------------------
# -------------- 35 V Vbias  ---------------
# ------------------------------------------

#EDIT next two lines: 470-515 for dark, 520-570 for LED
# --------------------------------------------------------------
vbias = 35; sipmN = '4'; source = 'led5'; #minR = 470; maxR = 515; 
peak = 'single'; connection = 'parallel'; header = 0; bins = 150; bins1 = 120;
# --------------------------------------------------------------
if source == 'dark': 
    minR = 470; maxR = 515;
else: 
    #minR = 560; maxR = 610;
    minR = 525; maxR = 570;  
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Set path to data  
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

#-------- Waveform plots ------------- 
#sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, header, vbias); # HEADER!!!!!!!!
phs_all, pis = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, header); # HEADER!!!!!!!!

#-------- Amplitude spectrum --------- 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins, -0.25, 0., vbias, connection, sipmN);
# ---------------------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
#plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[(len(ind_valleys)-5)]; #--------EDIT--------- 
hist_upperBound = ind_valleys[(len(ind_valleys)-1)]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE5 = x[ind_peak[(len(ind_peak)-6)]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins1, xMin, 0, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
#plt.close()

# -------------- EDIT X-AXIS ---------------
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
x = x/mu1;
# ------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
#x = [y-y[ind_peak[(len(ind_peak)-1)]] for y in x]  
#mu1 = x[ind_peak[(len(ind_peak)-2)]];
#x = x/mu1;
# Single PE
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
A1 = y[ind_peak[(len(ind_peak)-2)]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[(len(ind_peak)-3)]]; #--------EDIT--------- 
A2 = y[ind_peak[(len(ind_peak)-3)]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(vbias)+' V');
plt.legend(); plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[90:104],y[90:104],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[66:80],y[66:80],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2))), 

#elif peak == 'double': 
    #fit1 = fits.dualGaussFit(x[80:100],y[80:100],x, y,[mu1,sig,A1,(mu1-0.01),sig,(A1/3)]); #mu1 = delt[0]; sig1 = delt[1];
    #fit2 = fits.dualGaussFit(x[62:69],y[54:79],x, y,[mu2,sig,A2,(mu2-0.01),sig,(A2*2/3)]); #mu2 = delt[0]; sig2 = delt[1];
    #print('resolution: '+str(sig1/abs(mu1-mu2))), 

#plt.figure(); plt.plot(y)


# ------------------------------------------
# -------------- 35 V Vbias  ---------------
# ------------------------------------------

#EDIT next two lines: 470-515 for dark, 520-570 for LED
# --------------------------------------------------------------
vbias = 35; sipmN = '4'; source = 'dark'; #minR = 470; maxR = 515; 
peak = 'single'; connection = 'parallel'; header = 1; bins = 150; bins1 = 120;
# --------------------------------------------------------------
if source == 'dark': 
    minR = 470; maxR = 515;
else: 
    #minR = 560; maxR = 610;
    minR = 525; maxR = 570;  
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Set path to data  
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

#-------- Waveform plots ------------- 
#sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, header, vbias); # HEADER!!!!!!!!
phs_all, pis = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, header); # HEADER!!!!!!!!

#-------- Amplitude spectrum --------- 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins, -0.25, 0., vbias, connection, sipmN);
# ---------------------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
#plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[(len(ind_valleys)-6)]; #--------EDIT--------- 
hist_upperBound = ind_valleys[(len(ind_valleys)-1)]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE5 = x[ind_peak[(len(ind_peak)-6)]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins1, xMin, 0, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
#plt.close()

# -------------- EDIT X-AXIS ---------------
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
x = x/mu1;
# ------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
#x = [y-y[ind_peak[(len(ind_peak)-1)]] for y in x]  
#mu1 = x[ind_peak[(len(ind_peak)-2)]];
#x = x/mu1;
# Single PE
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
A1 = y[ind_peak[(len(ind_peak)-2)]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[(len(ind_peak)-3)]]; #--------EDIT--------- 
A2 = y[ind_peak[(len(ind_peak)-3)]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(vbias)+' V');
plt.legend(); plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[90:104],y[90:104],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[66:80],y[66:80],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2))), 

#elif peak == 'double': 
    #fit1 = fits.dualGaussFit(x[80:100],y[80:100],x, y,[mu1,sig,A1,(mu1-0.01),sig,(A1/3)]); #mu1 = delt[0]; sig1 = delt[1];
    #fit2 = fits.dualGaussFit(x[62:69],y[54:79],x, y,[mu2,sig,A2,(mu2-0.01),sig,(A2*2/3)]); #mu2 = delt[0]; sig2 = delt[1];
    #print('resolution: '+str(sig1/abs(mu1-mu2))), 

#plt.figure(); plt.plot(y)
"""

# ------------------------------------------
# -------------- 34 V Vbias  ---------------
# ------------------------------------------

#EDIT next two lines: 470-515 for dark, 520-570 for LED
# --------------------------------------------------------------
vbias = 35; sipmN = '4'; source = 'led2'; #minR = 470; maxR = 515; 
peak = 'single'; connection = 'parallel'; header = 0; bins = 150; bins1 = 120;
# --------------------------------------------------------------
if source == 'dark': 
    minR = 470; maxR = 515;
else: 
    #minR = 560; maxR = 610;
    minR = 525; maxR = 570;  
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Set path to data  
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

#-------- Waveform plots ------------- 
#sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, header, vbias); # HEADER!!!!!!!!
phs_all, pis = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, header); # HEADER!!!!!!!!

#-------- Amplitude spectrum --------- 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins, -0.25, 0., vbias, connection, sipmN);
# ---------------------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
#plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[(len(ind_valleys)-6)]; #--------EDIT--------- 
hist_upperBound = ind_valleys[(len(ind_valleys)-1)]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE5 = x[ind_peak[(len(ind_peak)-6)]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins1, xMin, 0, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
#plt.close()

# -------------- EDIT X-AXIS ---------------
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
x = x/mu1;
# ------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
#x = [y-y[ind_peak[(len(ind_peak)-1)]] for y in x]  
#mu1 = x[ind_peak[(len(ind_peak)-2)]];
#x = x/mu1;
# Single PE
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
A1 = y[ind_peak[(len(ind_peak)-2)]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[(len(ind_peak)-3)]]; #--------EDIT--------- 
A2 = y[ind_peak[(len(ind_peak)-3)]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(vbias)+' V');
plt.legend(); plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[88:106],y[88:106],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[67:82],y[67:82],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2))), 

#elif peak == 'double': 
    #fit1 = fits.dualGaussFit(x[80:100],y[80:100],x, y,[mu1,sig,A1,(mu1-0.01),sig,(A1/3)]); #mu1 = delt[0]; sig1 = delt[1];
    #fit2 = fits.dualGaussFit(x[62:69],y[54:79],x, y,[mu2,sig,A2,(mu2-0.01),sig,(A2*2/3)]); #mu2 = delt[0]; sig2 = delt[1];
    #print('resolution: '+str(sig1/abs(mu1-mu2))), 

#plt.figure(); plt.plot(y)
"""

# ------------------------------------------
# -------------- 34 V Vbias  ---------------
# ------------------------------------------

#EDIT next two lines: 470-515 for dark, 520-570 for LED
# --------------------------------------------------------------
vbias = 34; sipmN = '4'; source = 'dark'; #minR = 470; maxR = 515; 
peak = 'single'; connection = 'parallel'; header = 1; bins = 150; bins1 = 120;
# --------------------------------------------------------------
if source == 'dark': 
    minR = 470; maxR = 515;
else: 
    #minR = 560; maxR = 610;
    minR = 525; maxR = 570;  
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Set path to data  
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 

#-------- Waveform plots ------------- 
#sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, header, vbias); # HEADER!!!!!!!!
phs_all, pis = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, header); # HEADER!!!!!!!!

#-------- Amplitude spectrum --------- 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins, -0.25, 0., vbias, connection, sipmN);
# ---------------------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 15) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
#plt.close()

# Find peaks and valleys
hist_lowerBound = ind_valleys[(len(ind_valleys)-6)]; #--------EDIT--------- 
hist_upperBound = ind_valleys[(len(ind_valleys)-1)]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE5 = x[ind_peak[(len(ind_peak)-6)]]; 
# Replot with new range: first 4 peaks only 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, bins1, xMin, 0, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])

# Peaks and valleys for adjusted histo
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
#plt.close()

# -------------- EDIT X-AXIS ---------------
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
x = x/mu1;
# ------------------------------------------
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
#x = [y-y[ind_peak[(len(ind_peak)-1)]] for y in x]  
#mu1 = x[ind_peak[(len(ind_peak)-2)]];
#x = x/mu1;
# Single PE
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
A1 = y[ind_peak[(len(ind_peak)-2)]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[(len(ind_peak)-3)]]; #--------EDIT--------- 
A2 = y[ind_peak[(len(ind_peak)-3)]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.01; 
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ---------- MAKE PRETTY PLOT --------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(vbias)+' V');
plt.legend(); plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.show()

if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[88:106],y[88:106],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[67:82],y[67:82],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('RESOLUTION: '+str(sig1/abs(mu1-mu2))), 

#elif peak == 'double': 
    #fit1 = fits.dualGaussFit(x[80:100],y[80:100],x, y,[mu1,sig,A1,(mu1-0.01),sig,(A1/3)]); #mu1 = delt[0]; sig1 = delt[1];
    #fit2 = fits.dualGaussFit(x[62:69],y[54:79],x, y,[mu2,sig,A2,(mu2-0.01),sig,(A2*2/3)]); #mu2 = delt[0]; sig2 = delt[1];
    #print('resolution: '+str(sig1/abs(mu1-mu2))), 

#plt.figure(); plt.plot(y)
"""