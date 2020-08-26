"""
Full analysis script for 2019 SiPM paper 
Makes resolution fits for parallel connection 
"""
from matplotlib.ticker import FormatStrFormatter
import sipm_dataAnalysis
import detect_peaks
import matplotlib.pyplot as plt
import numpy as np
import pylab
import fits
 
# ------------------------------------------
# -------- ALWAYS SET PATH FIRST  ----------
# ------------------------------------------
# cd Documents/nEXO/sipm_analysis/code


# ------------------------------------------
# ----- Read out data and make plots  ------
# ------------------------------------------

#EDIT next two lines: 470-515 for dark, 520-570 for LED
# --------------------------------------------------------------
vbias = 65; sipmN = '4'; source = 'dark'; #minR = 470; maxR = 515; 
peak = 'single'; connection = 'series'; 
# --------------------------------------------------------------
if source == 'dark': 
    minR = 470; maxR = 515;
    header = 1;
else: 
    minR = 560; maxR = 610;
    header = 0; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
# Set path to data  
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
# Uncomment next line to get waveform plots 
#sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 10e6, header, vbias); # HEADER!!!!!!!!
phs_all, pis = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, header); # HEADER!!!!!!!!

# Integral spectrum -not good, commenting out for now
"""
plt.figure(); 
plt.hist(pis, 200, histtype = 'step', label = sipmN+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(vbias)+' V');
plt.legend();
#plt.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.xlabel('Pulse Integral (V)')
plt.ylabel('Events/bin')
    #plt.title('3 MHz BW filter, Vbias = '+str(Vbias)+' V')
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.show()
"""

# AMPLITUDE SPECTRUM #1: first PHS to find ranges 
#--------EDIT--------- 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, 100, -0.3, 0., vbias, connection, sipmN);
# Find peaks and valleys
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 10) # peaks in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) # valleys in PHS
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
print('Valleys: ', x[ind_valleys])
# Upper and lower bounds for Spectrum plots
hist_lowerBound = ind_valleys[(len(ind_valleys)-6)]; #--------EDIT--------- 
hist_upperBound = ind_valleys[(len(ind_valleys)-1)]
xMin = x[hist_lowerBound]; xMax = x[hist_upperBound];
PE5 = x[ind_peak[(len(ind_peak)-6)]]; 

# AMPLITUDE SPECTRUM #2: second PHS to find peaks for analysis 
x, y = sipm_dataAnalysis.get_PHSspectrum(phs_all, 100, xMin, xMax, vbias, connection, sipmN);
pylab.xlim([xMax,xMin])
# Edit mph and mpd IF necessary 
# peaks in PHS
ind_peak = detect_peaks.detect_peaks(y, mph = 20, mpd = 5) # CHANGE NEXT TWO WITH SAME MPD
# valleys in PHS
ind_valleys = detect_peaks.detect_peaks(y, mpd = 10, valley=True) 
plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); 
print('Valleys: ', x[ind_valleys])

"""
# ------------------------------------------
# -- ADJUST X-AXIS and make fit estimates --
# ------------------------------------------
mu0 = x[ind_peak[(len(ind_peak)-1)]];  # -1 or -2 depending on whether pedestal is detected
x = x/mu0;
# ------------------------------------------
# PEDESTAL
ind_peak2 = detect_peaks.detect_peaks(y, mph = 10, mpd = 5) # peaks in PHS
pedestal = x[ind_peak2[(len(ind_peak2)-1)]]
#pedestal = -0.02
 #  FOR UNDISTINGUISHABLE PEDESTAL SET MANUALLY 
x = x - pedestal # ADD +1 for indistinguishable pedestal   
# 5th PE for scaling 
PE4 = x[ind_peak2[(len(ind_peak2)-5)]]; # should match for indistinguishable pedestal
scale = 4/PE4
x = x*scale 
# ------------------------------------------
# Find fit estimates -ADJUST MPD APPROPRIATELY 
ind_peak = detect_peaks.detect_peaks(y, mph = 10, mpd = 5) # peaks in PHS
# Single PE
mu1 = x[ind_peak[(len(ind_peak)-2)]]; #--------EDIT--------- 
A1 = y[ind_peak[(len(ind_peak)-2)]];  #--------EDIT--------- 
# Second PE
mu2 = x[ind_peak[(len(ind_peak)-3)]]; #--------EDIT--------- 
A2 = y[ind_peak[(len(ind_peak)-3)]];  #--------EDIT--------- 
# Expected fit values 
sig = 0.05;
expected1 = [mu1,sig,A1]; 
expected2 = [mu2,sig,A2];

# ------------------------------------------
# ---------- MAKE PRETTY PLOT --------------
# ------------------------------------------
fig, ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.rcParams.update({'font.size': 12})
plt.plot(x[0:len(y)],y,'k.', label = sipmN+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(vbias)+' V');
plt.legend(); plt.xlabel('PE'); plt.ylabel('Events/bin'); plt.show()
#plt.plot(xPE[0:len(y)],y,'k.', label = sipmN+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(vbias)+' V');
#plt.legend(); plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin'); plt.show()
#pylab.xlim([xMax,xMin])
# ------------- FITS -------------
if peak == 'single': 
    delta = []; 
    # ------ EDIT FIT RANGES ------ speMin:speMax  pe2Min:pe2Max
    delt, sig = fits.singleFitter(x[84:93],y[84:93],x,y, expected1); mu1 = delt[0]; sig1 = delt[1];#delta.append(delt[0])
    delt, sig = fits.singleFitter(x[69:75],y[69:75],x,y, expected2); mu2 = delt[0]; sig2 = delt[1];#delta.append(delt[0])
    print('resolution: '+str(sig1/abs(mu1-mu2))), 
    

elif peak == 'double': 
    fit1 = fits.dualGaussFit(x[80:100],y[80:100],x, y,[mu1,sig,A1,(mu1-0.01),sig,(A1/3)]); #mu1 = delt[0]; sig1 = delt[1];
    fit2 = fits.dualGaussFit(x[68:74],y[68:74],x, y,[mu2,sig,A2,(mu2-0.01),sig,(A2*2/3)]); #mu2 = delt[0]; sig2 = delt[1];
    #print('resolution: '+str(sig1/abs(mu1-mu2))), 

plt.figure(); plt.plot(y)


#------------------------------------------------------------------------------
resolutions =[] ; nsipms = [] 
m,b = np.polyfit(nsipms,resolutions,1)
yList = []
for i in range(len(nsipms)):
    yFit = m*nsipms[i]+b
    yList.append(yFit)
plt.figure(); plt.scatter(nsipms, resolutions); plt.plot(nsipms, yList, linestyle='solid')
"""

"""
# DARK NOISE 
vbias = 62; minR = 470; maxR = 515; 
pathToData = 'sipm_data/4sipm_testing/series_test/wave0_'+str(vbias)+'V_dark.dat'; 
sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, 1, vbias); 
phs, phs2, phs3 = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, 1);  
x, y = sipm_dataAnalysis.get_PHSspectrum(phs, 200, -0.25, 0., vbias);

# LED 
vbias = 66; minR = 560; maxR = 610; 
pathToData = 'sipm_data/6sipm_testing/series_test/wave0_'+str(vbias)+'V_led.dat'; 
sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, 1, vbias); 
phs, phs2, phs3 = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, 1);  
x, y = sipm_dataAnalysis.get_PHSspectrum(phs, 200, -0.25, 0., vbias);
"""