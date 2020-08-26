"""
Full analysis script for 2019 SiPM paper. Makes ENC fits for 2-SiPMs 
(2s) for t-rise dependence study, using the triangular filter. Uses sample 
biased at Vbias = -68 V (Vover ~ 4.5 V)
"""
# Set path to SiPM methods
# insert at 1, 0 is the script path (or '' in REPL)
import sys
sys.path.insert(1, '/Users/johny/Documents/nEXO/sipm_analysis/code/methods')

# Import all other needed libraries 
from signal_processing import get_amplitudes
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import detect_peaks
import fits
import time

# Start clock. Used during optimization 
start_time = time.time()

# Local path: 
# cd Documents/nEXO/sipm_analysis/code/scripts

## ################## ##
##     MAIN BODY      ##
## ################## ##

if __name__ == "__main__":
    
    # Read data file from Mietek: 
    #df_parallel_gauss = pd.read_csv('simulations_mietek/enc_semigauss_shaper_parallel.csv')
    df_series_triang = pd.read_csv('../noise-model/enc-triangFilter-noiseModel.csv')
    tpeak_series_triang = df_series_triang['peaking time (ns)'].tolist()
    enc_series_triang_6sipm = df_series_triang['ENC-2s'].tolist()
    
    # Initiate lists for observables:
    resolutions_2s, errors, Tlist_2s  = [], [], []
    # EXPECTED VALUES 
    mu1 = 1; A1 = 0.3; sig = 0.1; mu2 = 2; A2 = 0.02;  
    expected1 = [mu1,sig,A1]; expected2 = [mu2,sig,A2];
    
    # Set global variables: 
    Npulses = 20000; vbias = 68; connection = 'series'; peak = 'single';
    peMin_unfiltered = 0.02; 
    
    # 6-sipms 2s3p  
    sipmN = '2'; source = 'dark'; header = 1;
    pathToData = '../../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    plt.ion()
    
    
    # ----------------------------- 25 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 25; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 1), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 19, 28; 
    xR, yR = x[xMin:xMax], y[xMin:xMax]
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    #plt.close() 
    
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 50 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 50; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 2), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 17, 25; 
    xR, yR = x[xMin:xMax], y[xMin:xMax]
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 75 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 75; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 2.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 18, 27; 
    xR, yR = x[xMin:xMax], y[xMin:xMax]
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 100 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 100; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 18, 26; 
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 100 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 125; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 20, 29; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    plt.close() 
    
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 150 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 150; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 32; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    plt.figure(); plt.plot(y)
    
    # ----------------------------- 175 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 175; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 20, 33; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 200 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 200; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 35; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 225 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 225; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 22, 36; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 225 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 250; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 22, 37; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 275 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 275; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 22, 38; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 300; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[1]];
    x = x/mu1;
    
    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 22, 39; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 350; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 36; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 400; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 36; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)

    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 450; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 36; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 500; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 37; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)

    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 600; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 38; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 700; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 40; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)

    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 800; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 21, 41; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)
    
    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 900; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 18, 43; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)

    # ----------------------------- 300 ns -----------------------------------------
    # Set peaking time:  
    peMin = 0.015; peMax = 0.028; peakingTime = 1000; Tlist_2s.append(peakingTime)
    
    # Compute rise time and get amplitudes 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, t_peak = peakingTime, filtering = 'triangular'); 
    
    # Plot PE spectrum:
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3.5), density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    #y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', range = (0, 3), label = 'V$_{bias}$ = -'+str(vbias)+' V, $t_p$ = '+str(peakingTime)+' ns');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    
    # Find peaks and valleys #1 
    ind_peak = detect_peaks.detect_peaks(y, mpd = 10) # peaks in PHS
    ind_valleys = detect_peaks.detect_peaks(y, mpd = 15, valley=True) # valleys in PHS
    plt.plot(x[ind_peak],y[ind_peak], 'o', color='black'); #TEST
    plt.plot(x[ind_valleys],y[ind_valleys], 'o', color='red'); #TEST
    #plt.close()
    
    # Normalize by dividing by first PE peak value to convert mV -> PE 
    mu1 = x[ind_peak[2]];
    x = x/mu1;

    # Make initial fit: 
    delta = []; 
    # RANGES: 
    xMin, xMax = 18, 44; #!!!
    delt, sig = fits.singleFitter(x[xMin:xMax], y[xMin:xMax], expected1, plotting = 0); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); #delta.append(delt[0]) 
    # Correct axis and make secondary fit
    xAxis = np.linspace(0,2,50)
    gaussianFit = fits.gauss(x, mu1, sig1, A1)
    gaussianFit = list(gaussianFit)
    index = gaussianFit.index(max(gaussianFit))
    x = x/(x[index])
    gaussianFit = fits.gauss(xAxis, mu1, sig1, A1)
    #plt.figure(); plt.plot(y)
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    plt.plot(x[0:len(y)],y,'k.', label = 'V$_{over}$ = '+str((vbias-59)/2.)+' V, '+'Gauss filter, $t_{peak} = $: '+str(peakingTime));
    plt.xlabel('PE'); plt.ylabel('Normalized counts'); plt.title(sipmN+' SiPMs: 2s')
    plt.legend(loc='upper right');  plt.show()
    #plt.plot(xAxis, gaussianFit, color='red')
    
    delta = []; 
    delt, error = fits.singleFitter(x[xMin:xMax],y[xMin:xMax], expected1); 
    mu1 = delt[0]; sig1 = round(delt[1], 3); A1 = delt[2]; enc_error = round(error[1], 3)
    print('ENC for Tp = '+str(peakingTime)+' ns: '+str(sig1)+' +/-'+str(enc_error)+' s.p.e.')
    errors.append(enc_error) 
    
    resolutions_2s.append(abs(sig1));
    #resolutions_2s.append(abs(sig1)/abs(mu1-mu2)); 
    plt.xlim(0, 4);
    
    #plt.figure(); plt.plot(y)

    
    plt.close('all')
    
    # ------------------------------------------------------------------------------
    # ---------------------------- ENC plot ---------------------------------
    # ------------------------------------------------------------------------------
    
    # Make ENC vs Tpeak plot: 
    plt.figure(); 
    Tlist_2s = [x*1e-9 for x in Tlist_2s] 
    plt.errorbar(Tlist_2s, resolutions_2s, yerr=errors, fmt='o', capsize = 2, capthick=2, label = 'Data')
    #plt.plot(Tlist_2s, resolutions_2s, 'b.', label = 'Data')
    plt.plot(tpeak_series_triang, enc_series_triang_6sipm, 'r', label = 'Model')
    plt.title('Triangular filter, 2 SiPMs: 2s, $V_{over}$ = 4.5 V') ; 
    plt.xlabel('$t_{peak}$ (s)'); plt.ylabel('enc ( s.p.e.)');
    plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()
    plt.savefig('../figures/enc-tp/triangular/enc-triangFilter-2s.png')
    
    # Make pandas dataframes: 
    df_freq_2s = pd.DataFrame({'peaking time (ns)': Tlist_2s,
                    'ENC(s.p.e.)': resolutions_2s,
                    'error(s.p.e.)': errors})
                    
    # Save as .csv file: 
    df_freq_2s.to_csv('../data-output/enc-tp/triangular/enc-triangFilter-2s.csv')
    
    print("--- %s seconds ---" % (time.time() - start_time))
    input('Press [ENTER] to finish.')