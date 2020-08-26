"""
A set of methods to analyze the nuisance parameters from FBK SiPMs using 
data from a CAEN digitizer. 
Methods for dark noise, afterpulsing, and computing the charge in mV*ns.
"""

from numpy import fromfile, dtype, linspace, mean
import matplotlib.pyplot as plt
import numpy as np
import wavedumpReader
import detect_peaks
import BWfilter
import math

# ------------------------------------------------------------------------------
# ---------------------------- ALWAYS SET PATH FIRST  --------------------------
# ------------------------------------------------------------------------------
# cd Documents/nEXO/sipm_analysis/code

def getTimeArray(fileName, PEmin, PEmax, Nsipms, Vbias, cutoffFreq = 10e6):
    """
    A routine to extract the time stamp array from headers in CAEN WaveDump files.
    This is then used to compute the dark noise. 

    Parameters: 
    ------------
    - fileName: input file 
    - singlePEcut: cut out the pedestal and only pick > 1 PE events, use absolute val
    - Nsipms: number of SiPMs
    - Vbias: biasing voltage

    Returns: 
    ------------
    timeList: array of time stamps for each trigger

    Example use: 
    ------------
    See getDarkNoise method

    Notes: 
    ------------
    This method only takes in data files ***with header***
    """
    dataFile = wavedumpReader.DataFile(fileName)
    dataFile.file.seek(0) 
    Npulses, timeTagRollover, i, oldTimeTag = 0, 0, 0, 0.  
    NpulseList, timeList, timeDiff = [], [], []
    while True:#for i in range(N):
        header = fromfile(dataFile.file, dtype='I', count=6)
        if len(header) != 6:
            break
        eventSize = (header[0] - 24) // 2
        trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize)#count=1012
        trace = trace - mean(trace[0:100]);
        trace = trace*2/(2**12);
        T = linspace(0,4*len(trace),len(trace))
        fs = 256e6;
        filteredSig = BWfilter.butter_lowpass_filter(trace,cutoffFreq,fs)#
        ind_valleys = detect_peaks.detect_peaks(filteredSig, mph = PEmin, valley=True)
        #---------------------
        triggerTimeTag = header[5]
        if triggerTimeTag < oldTimeTag:
            timeTagRollover += 1
            oldTimeTag = float(header[5])
        else:
            oldTimeTag = float(header[5])
        # Correcting triggerTimeTag for rollover
        triggerTimeTag += timeTagRollover*(2**31)
        # Convert from ticks to ns since the beginning of the file
        triggerTime = triggerTimeTag * 8 
        if abs(min(filteredSig)) > abs(PEmin): #and abs(PEmax) > abs(min(filteredSig)):
            Npulses += 1 
            timeList.append(triggerTime) 
            NpulseList.append(Npulses)
            #i += 1; 
            if i >= 1: 
                delta_t = triggerTime - timeList[i-1]
                timeDiff.append(delta_t)
            i += 1
            if len(ind_valleys) > 1:
                for j in range(len(ind_valleys)):
                    delta_t = T[ind_valleys[j]] - T[ind_valleys[0]]
                    #apTime = triggerTime + timeDiff 
                    timeDiff.append(delta_t)
    # Sanity check for time list, should increase mostly linearly
    #plt.figure(); plt.plot(timeList)
    return timeDiff
    
def get_charge(fileName, PEmin, PEmax, cutoffFreq = 10e6):
    """
    A method to extract the average single PE charge from 
    CAEN wavedump files. 

    Parameters: 
    ------------
    - fileName: input file 
    - PEmin: lower limit for SPE 
    - PEmax: lower limit for SPE

    Returns: 
    ------------
    pulse: average of 1000 pulses, rise time shown in plot

    Example use: 
    ------------
    from sipm_signalProcessing import getAvgPulse
    peMin = 0.021; peMax = 0.06; 
    vbias = 71; sipmN = '4'; source = 'dark'; connection = 'series'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
    charge = get_charge(pathToData, peMin, peMax, sipmN, vbias)

    Notes: 
    ------------
    This method only takes in data files with header
    """
    dataFile = wavedumpReader.DataFile(fileName)
    dataFile.file.seek(0) 
    i = 0; pulses = [];  
    while True:
        header = fromfile(dataFile.file, dtype='I', count=6)
        if len(header) != 6:
            break
        eventSize = (header[0] - 24) // 2
        trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize)
        trace = trace - mean(trace[0:100]);
        trace = trace*2/(2**12);
        T = linspace(0,4*len(trace),len(trace))
        fs = 256e6;
        filteredSig = BWfilter.butter_lowpass_filter(trace,cutoffFreq,fs)
        ind_valleys = detect_peaks.detect_peaks(filteredSig, mph = PEmin, valley=True)
        #---------------------
        if abs(PEmax) > abs(min(filteredSig)) and abs(min(filteredSig)) > abs(PEmin):
            if len(ind_valleys) == 1:
                pulses.append(filteredSig);
                i += 1;
                if i == 1000: 
                    break
    pulse_neg = 1000*(np.mean(pulses, axis=0))
    #plt.figure(); plt.plot(T, pulse_neg, label = 'negative'); plt.show(); plt.legend();
    pulse_pos = [-x for x in pulse_neg]
    plt.figure(); plt.plot(T, pulse_pos, label = 'positive'); plt.show(); plt.legend();
    charge =  np.trapz(pulse_pos, x=T);
    return charge

def getDarkNoise_test(timeArray, Nsipms, Vbias, binNum = 200, MIN = 1, MAX = 1e9): 
    """
    A method to compute the dark noise rate. Takes inputs from getTimeArray, 
    and produces distribution plot for time difference between secondary 
    and primary pulses. From this, the dark noise is found by fitting the 
    linear tail of the distribution.   

    Parameters: 
    ------------
    - timeList: time stamps array, output from getArrays 
    - Npulses: raw number of dark pulses, output from getArrays
    - Vbias: Biasing voltage, input set directly 

    Returns: 
    ------------
    dnRateAvg: avg dark noise rate for a single data file 
    errorBar: error bar corresponding to the measurement 

    Example use: 
    ------------
    from nuissanceParams_methods import getDarkNoise, getTimeArray 
    peMin = 0.018; peMax = 0.054; 
    vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
    time = getTimeArray(pathToData, peMin, peMax, sipmN, vbias)
    getDarkNoise_test(time, sipmN, vbias)
    #plt.yscale('log')
    """
    #timeList = []
    #for i in range(len(timeArray)-1): 
        #timeDiff = (timeArray[i+1]-timeArray[i])
        #timeList.append(timeDiff) 
    # Normalized histogram: 
    plt.figure();
    y,x,_= plt.hist(timeArray, bins=np.logspace(np.log10(MIN),np.log10(MAX), binNum), histtype='step', normed=True, label=str(Nsipms)+' SiPMs biased at '+str(Vbias)+' V')
    time_ax = []
    for i in range(len(x)-1):
        t = (x[i+1]+x[i])/2
        time_ax.append(t)
    lamb_list, rho_list = [], []
    lamb_list.append(0.); 
    for i in range(len(y)): 
        beta = sum(lamb_list) 
        lambd = -np.log(1-y[i]*(timeArray[i+1]-timeArray[i])/math.exp(-beta))
        #lambd = -np.log(1-y[i]/math.exp(-beta))
        lamb_list.append(lambd)
        delta_t = (x[i+1]-x[i])
        rho = (lambd/delta_t)*10**9
        rho_list.append(rho)
    #binSize = (max(x)-min(x))/binNum;
    plt.yscale('log'); plt.xlabel('time difference (ns)'); plt.ylabel('Events/bin'); 
    plt.legend(); plt.gca().set_xscale("log"); #plt.xscale('log');
    #plt.close()
    rho_new = [x/(float(Nsipms)*100.) for x in rho_list] 
    plt.figure(); plt.plot(time_ax,rho_new,'.', label=' SiPMs biased at '+str(Vbias)+' V'); plt.yscale('log'); plt.xscale('log')
    plt.xlabel('time difference (ns)'); plt.ylabel('Hz/mm$^2$'); plt.legend()
# Test

def getDarkNoise(timeArray, Nsipms, Vbias, binNum = 200, MIN = 1, MAX = 1e9): 
    """
    A method to compute the dark noise rate. Takes inputs from getTimeArray, 
    and produces distribution plot for time difference between secondary 
    and primary pulses. From this, the dark noise is found by fitting the 
    linear tail of the distribution.   

    Parameters: 
    ------------
    - timeList: time stamps array, output from getArrays 
    - Npulses: raw number of dark pulses, output from getArrays
    - Vbias: Biasing voltage, input set directly 
    
    Returns: 
    ------------
    dnRateAvg: avg dark noise rate for a single data file 
    errorBar: error bar corresponding to the measurement 
    
    Example use: 
    ------------
    from nuissanceParams_methods import getDarkNoise, getTimeArray 
    peMin = 0.018; peMax = 0.054; 
    vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
    time = getTimeArray(pathToData, peMin, peMax, sipmN, vbias)
    getDarkNoise(time, sipmN, vbias)
    #plt.yscale('log')
    """
    timeList = []
    for i in range(len(timeArray)-1): 
        timeDiff = (timeArray[i+1]-timeArray[i])
        timeList.append(timeDiff)
    # Make histogram 
    plt.figure(); 
    #binNum = 200; xmin = 1; xmax = 1e8; 
    #binSize = (xmax-xmin)/binNum; 
    #y,x,_= plt.hist(timeList, range = (xmin,xmax), bins = binNum, histtype='step', label=str(Nsipms)+' SiPMs biased at '+str(Vbias)+' V'); 
    y,x,_= plt.hist(timeList, bins=np.logspace(np.log10(MIN),np.log10(MAX), binNum), histtype='step', normed=True, label=str(Nsipms)+' SiPMs biased at '+str(Vbias)+' V')
    plt.yscale('log'); plt.xlabel('time difference (ns)'); plt.ylabel('Events/bin'); plt.legend(); plt.gca().set_xscale("log"); #plt.xscale('log');
    plt.close()
    # Rescale to obtain y-axis as Hz/mm^2 and plot
    #y = (y*10**9)/(binSize*float(Nsipms)*100.)
    y = (y*10**9)/(float(Nsipms)*100.)
    plt.figure(); plt.plot(x[0:len(y)],y,'.', label=' SiPMs biased at '+str(Vbias)+' V'); plt.yscale('log'); plt.xscale('log')
    plt.xlabel('time difference (ns)'); plt.ylabel('Hz/mm$^2$'); plt.legend()
# Test 

def afterpulsing_test(inputFile, PEmin, PEmax, Header, Vbias, cutoffFreq = 10e6):    
    """
    A method to plot a sample of N pulses from a CAEN wavedump file. 

    Parameters:
    ----------- 
    inputFile: string, input file-name
    N: int, number of events to be read 
    Header: 1 if input file has header, 0 otherwise 
    Vbias: biasing voltage 
    cutoffFreq: cutoff frequency in Hertz, 10 MHz is best for sipm paper

    Example use:
    -----------
    from sipm_signalProcessing import afterpulsing_test 
    Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
    if source == 'dark': 
        minR = 470; maxR = 515;
        header = 1;
    else: 
        minR = 560; maxR = 610;
        header = 0; 
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias);
    nbins = 200; lower = 0.; upper = 80; 
    plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
    """
    dataFile = wavedumpReader.DataFile(inputFile)
    dataFile.file.seek(0)
    plt.rcParams.update({'font.size': 12})
    AP = []; i = 0; amplitudes_pe = []; amplitudes = []
    # Calculate average charge from 1000 waveforms w/no afterpulsing
    Q0 = get_charge(inputFile, PEmin, PEmax)
    # Scale up PE limits to mV
    PEmin = 1000*PEmin; PEmax = 1000*PEmax; 
    # Dark noise: 
    if Header == 1: 
    #while True: 
        minR = 470; maxR = int(2200/4); full = int(minR+1000/4); 
        #for i in range(N):
        while True: 
            header = fromfile(dataFile.file, dtype='I', count=6)
            if len(header) != 6:
                break
            i += 1
            eventSize = (header[0] - 24) // 2
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            fs = 256e6
            filteredSig = (BWfilter.butter_lowpass_filter(trace,cutoffFreq,fs))
            pulse = 1000*([-x for x in filteredSig])
            timeAx = linspace(0,4*len(trace),len(trace))
            #if i == 1:
                #plt.figure(); plt.plot(timeAx, pulse)
            #amplitude = max(pulse[minR:maxR]); amplitude_list.append(amplitude)
            if abs(PEmax) > abs(max(pulse)) and abs(max(pulse)) > abs(PEmin):
                #Q0 = np.trapz(pulse[minR:maxR], x=timeAx[minR:maxR]); #integral_list.append(integral)
                Qt = np.trapz(pulse[minR:full], x=timeAx[minR:full]);
                A = Qt/Q0 -1
                AP.append(A)
                amplitudes_pe.append(max(pulse))
                # Test: 
            if i == 1: 
                plt.figure(); plt.plot(timeAx, pulse) 
                #plt.axvline(x=timeAx[minR], color='r',linestyle='--', label = str(timeAx[minR]))
                #plt.axvline(x=timeAx[maxR], color='r',linestyle='--', label = str(timeAx[minR]))
                #plt.axvline(x=timeAx[full], color='r',linestyle='--', label = str(timeAx[full]))
            amplitudes.append(max(pulse))
    # LED trigger:  
    else: 
        minR = 560; maxR = 610;
        for i in range(N):
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=1024) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            fs = 256e6;
            filteredSig = 1000*(BWfilter.butter_lowpass_filter(trace,cutoffFreq,fs))
            timeAx = linspace(0,4*len(trace),len(trace))
    del dataFile;
    #plt.hist(amplitudes, 200, histtype = 'step', label = 'Full spectrum');
    plt.figure(); plt.hist(amplitudes_pe, 200, histtype = 'step', label = '1 PE');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.figure(); plt.hist(amplitudes, 200, histtype = 'step', label = 'Full spectrum');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    afterpulsing = np.mean(AP)
    print('Afterpulsing rate for '+str(Vbias)+'V is: ', afterpulsing)
    return afterpulsing
   
def getAfterpulsing(fileName, PEmin, PEmax, Nsipms, Vbias, cutoffFreq = 10e6): #Use the same dataFile defined above
    """
    A routine to extract the time stamp array from headers in CAEN WaveDump files.
    This is then used to compute the dark noise. 

    Parameters: 
    ------------
    - fileName: input file 
    - singlePEcut: cut out the pedestal and only pick > 1 PE events, use absolute val
    - Nsipms: number of SiPMs
    - Vbias: biasing voltage

    Returns: 
    ------------
    timeList: array of time stamps for each trigger

    Example use: 
    ------------
    peMin = 0.011; peMax = 0.03; 
    vbias = 65; sipmN = '4'; source = 'dark'; connection = 'series'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
    getAfterpulsing(pathToData, peMin, peMax, sipmN, vbias)

    Notes: 
    ------------
    This method only takes in data files ***with header***
    """
    dataFile = wavedumpReader.DataFile(fileName)
    dataFile.file.seek(0) 
    AP = []
    Q0 = get_charge(fileName, PEmin, PEmax)
    while True:
    #for i in range(20):
        header = fromfile(dataFile.file, dtype='I', count=6)
        if len(header) != 6:
            break
        eventSize = (header[0] - 24) // 2
        trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize)
        trace = trace - mean(trace[0:100]);
        trace = trace*2/(2**12);
        T = linspace(0,4*len(trace),len(trace))
        fs = 256e6;
        filteredSig = BWfilter.butter_lowpass_filter(trace,cutoffFreq,fs)
        ind_valleys = detect_peaks.detect_peaks(filteredSig, mph = PEmin, valley=True)
        # 1 P.E. cut 
        if abs(min(filteredSig)) > abs(PEmin) and abs(PEmax) > abs(min(filteredSig)):
            #i += 1; 
            Q0 = filteredSig[ind_valleys[0]]
            #Qt = sum(filteredSig[ind_valleys])
            if len(ind_valleys) == 1:
                A = 0.
                AP.append(A) 
            if len(ind_valleys) > 1:
                #for j 
                #if (T) < 1000
                #Qt = sum(filteredSig[ind_valleys])
                #A = Qt/Q0 -1
                #AP.append(A) 
                Qlist = []
                for j in range(len(ind_valleys)-1):
                    timeDiff = T[ind_valleys[j+1]] - T[ind_valleys[0]]
                    if timeDiff < 1000: 
                        Qlist.append(filteredSig[ind_valleys[j+1]])
                    #apTime = triggerTime + timeDiff 
                    #timeList.append(apTime)
                    # Sanity check/debugging
                    #plt.plot(T[ind_valleys[j+1]],filteredSig[ind_valleys[j+1]], 'o', color='red');
                Qt = sum(Qlist)+Q0
                A = Qt/Q0 -1
                AP.append(A)
            # Sanity check/debugging
            #if i < 20: 
                #plt.plot(T, filteredSig)
                #plt.plot(T[ind_valleys],filteredSig[ind_valleys], 'o', color='red');
    afterpulsing = mean(AP)
    print('Afterpulsing rate for '+str(Vbias)+'V is: ', afterpulsing)
        
def plotDN(vArray, dnArray, dnErrorbars): 
    """
    A method to plot the dark noise rate as a function of overvoltage
    
    Parameters: 
    ------------
    vArray: overvoltage array 
    dnArray: dark noise array. 
    dnErrorbars: error bar array
    All arrays must be the same size. 
    
    Returns: 
    ------------
    Plot of dark noise (y-axis) as a function of overvoltage (x-axis) with error bars
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 13})
    plt.figure(); plt.errorbar(vArray, dnArray, yerr = dnErrorbars, fmt='o', capsize = 2, capthick=2); plt.grid(True)
    plt.xlabel('Overvoltage (V)')
    plt.ylabel('Dark Noise Rate [Hz/mm$^2$]')
