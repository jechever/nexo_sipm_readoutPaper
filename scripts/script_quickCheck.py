"""
Quick analysis script for checking data for 2019 SiPM paper.  
Produces plot with 100 sample waveforms and pulse height spectrum. 
Need only to modify lines 15 and 22 according to biasing voltage of file,
range is already optimized for dark data files and LED files. 
Make sure to save files in same format and modify path to file accordingly, 
e.g. for a file with dark noise data at 68 V file name should be 
wave0_68V_dark.dat
"""

import sipm_dataAnalysis

 
# ------------------------------------------
# ----- Read out data and make plots  ------
# ------------------------------------------

# Set path in command line
# cd Documents/nEXO/sipm_analysis

# Set parameters
vbias = 66; nSipms = 6; connection = 'series'; 

# DARK NOISE 
minR = 470; maxR = 515; 
pathToData = '../sipm_data/'+str(nSipms)+'sipm_testing/series_test/wave0_'+str(vbias)+'V_dark.dat'; 
sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, 1, vbias); 
phs, phs2 = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, 1);  
x, y = sipm_dataAnalysis.get_PHSspectrum(phs, 200, -0.25, 0., vbias, connection, nSipms);

# LED 
minR = 560; maxR = 610; 
pathToData = 'sipm_data/'+str(nSipms)+'sipm_testing/series_test/wave0_'+str(vbias)+'V_led.dat'; 
sipm_dataAnalysis.get_filteredPlots(pathToData, 100, 0.03, 1, vbias); 
phs, phs2 = sipm_dataAnalysis.getPHS(pathToData, 20000, 0.03, 0.02, minR, maxR, 1);  
x, y = sipm_dataAnalysis.get_PHSspectrum(phs, 200, -0.25, 0., vbias, connection, nSipms);
