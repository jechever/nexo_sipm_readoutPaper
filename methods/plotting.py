from numpy import fromfile, dtype, linspace
from scipy import ndimage
import matplotlib.pyplot as plt
import numpy as np 
import wavedumpReader
import filters      

#-------------------------------------------------------------------------------
#--------------------------- Plotting methods ----------------------------------
#-------------------------------------------------------------------------------

def plot_peSpectrum(PHS, Nbins, minVal, maxVal, Vbias, connection, nSIPMS): 
    """
    A function to plot the pulse height spectrum (PHS) from a SiPM data set
    
    Parameters:
    ----------- 
    PHS: array of pulse amplitudes
    Nbins: number of bins 
    minVal: minimum value you want the histogram to have
    maxVal: maximum value you want the histogram to have
    Vbias: biasing voltage
    
    Returns: 
    Plot: pulse height spectrum (PHS) histogram 
    yf: y-values of PHS
    xf: x-values of PHS
          
    Example use:
    ------------ 
    from signal_processing import
    y33, x33 = plot_peSpectrum(PHS33, 100, -0.20, 0.0, 33);
    """
    #from matplotlib.ticker import FormatStrFormatter
    fig, ax = plt.subplots()
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 12})
    yf,xf,_= plt.hist(PHS, Nbins, histtype = 'step', range = (minVal,maxVal), label = nSIPMS+' SiPMs, '+connection+' connection'+', V$_{bias}$ = -'+str(Vbias)+' V');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()
    #plt.title('3 MHz BW filter, Vbias = '+str(Vbias)+' V')
    plt.show()
    return xf, yf
    
def plot_rawPulses(inputFile, Header, N = 100):    
    """
    A function to return the raw (unfiltered) pulses plot, takes 
    in binary files from CAEN. 
    
    Parameters:
    ------------ 
    
    Example use:
    ------------
    from signal_processing import plot_rawPulses
    vbias = 65; sipmN = '4'; source = 'dark'; 
    peak = 'single'; connection = 'series'; 
    if source == 'dark': 
        minR = 470; maxR = 515;
        header = 1;
    else: 
        minR = 560; maxR = 610;
        header = 0; 
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    plot_rawPulses(pathToData, header, vbias);     
    """
    plt.rcParams.update({'font.size': 12})
    dataFile = wavedumpReader.DataFile(inputFile)
    dataFile.file.seek(0)
    plt.figure()
    if Header == 1: 
        minR = 470; maxR = 515;
        for i in range(N):
            header = fromfile(dataFile.file, dtype='I', count=6)
            if len(header) != 6:
                break
            eventSize = (header[0] - 24) // 2
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
            #trace = trace - np.mean(trace[0:100]);
            #trace = 1000*(trace*2/(2**12));
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.plot(timeAx, trace);
    # LED trigger:  
    else: 
        minR = 560; maxR = 610;
        for i in range(N):
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=1024) 
            #trace = trace - np.mean(trace[0:100]);
            #trace = 1000*(trace*2/(2**12));
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.plot(timeAx, trace); 
    #plt.axvline(x=timeAx[minR], color='r',linestyle='--')
    #plt.axvline(x=timeAx[maxR], color='r',linestyle='--')
    del dataFile;
    plt.xlabel('Time (ns)'); plt.ylabel('a.u. (ADC)')
    plt.title('Raw signal')
    return trace

def plot_filteredPulses(inputFile, Header, Vbias, cutoffFreq = 10e6, N = 100):    
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
    from signal_processing import plot_filteredPulses
    vbias = 65; sipmN = '4'; source = 'dark'; 
    peak = 'single'; connection = 'series'; 
    if source == 'dark': 
        minR = 470; maxR = 515;
        header = 1;
    else: 
        minR = 560; maxR = 610;
        header = 0; 
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    plot_filteredPulses(pathToData, header, vbias); 
    """
    dataFile = wavedumpReader.DataFile(inputFile)
    dataFile.file.seek(0)
    plt.rcParams.update({'font.size': 12})
    plt.figure() 
    # Dark noise: 
    if Header == 1: 
        minR = 470; maxR = 515;
        for i in range(N):
            header = fromfile(dataFile.file, dtype='I', count=6)
            if len(header) != 6:
                break
            eventSize = (header[0] - 24) // 2
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            fs = 256e6
            filteredSig = 1000*(filters.butter_lowpass_filter(trace,cutoffFreq,fs))
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.plot(timeAx, filteredSig);
    # LED trigger:  
    else: 
        minR = 525; maxR = 610;
        for i in range(N):
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=1024) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            fs = 256e6;
            filteredSig = 1000*(filters.butter_lowpass_filter(trace,cutoffFreq,fs))
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.plot(timeAx, filteredSig); 
    plt.axvline(x=timeAx[minR], color='r',linestyle='--')
    plt.axvline(x=timeAx[maxR], color='r',linestyle='--')
    del dataFile;
    plt.xlabel('Time (ns)'); plt.ylabel('Voltage (mV)')
    plt.title('10 MHz BW filter, $V_{bias}$ = -'+str(Vbias)+' V')
    return filteredSig

def plot_filteredPulses_gauss(inputFile, Header, Vbias, t_peak = 70, N = 100):    
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
    from signal_processing import plot_filteredPulses
    vbias = 65; sipmN = '4'; source = 'dark'; 
    peak = 'single'; connection = 'series'; 
    if source == 'dark': 
        minR = 470; maxR = 515;
        header = 1;
    else: 
        minR = 560; maxR = 610;
        header = 0; 
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    plot_filteredPulses(pathToData, header, vbias); 
    """
    dataFile = wavedumpReader.DataFile(inputFile)
    dataFile.file.seek(0)
    plt.rcParams.update({'font.size': 12})
    plt.figure() 
    # Dark noise: 
    if Header == 1: 
        minR = 450; maxR = 515;
        sigma = t_peak/15
        for i in range(N):
            header = fromfile(dataFile.file, dtype='I', count=6)
            if len(header) != 6:
                break
            eventSize = (header[0] - 24) // 2
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            filteredSig = 1000*(ndimage.gaussian_filter1d(trace, sigma))
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.plot(timeAx, filteredSig);
    # LED trigger:  
    else: 
        minR = 500; maxR = 610;
        sigma = t_peak/15
        for i in range(N):
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=1024) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            filteredSig = 1000*(ndimage.gaussian_filter1d(trace, sigma))
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.plot(timeAx, filteredSig); 
    plt.axvline(x=timeAx[minR], color='r',linestyle='--')
    plt.axvline(x=timeAx[maxR], color='r',linestyle='--')
    del dataFile;
    plt.xlabel('Time (ns)'); plt.ylabel('Voltage (mV)')
    plt.title('Gaussian filter filter, $V_{bias}$ = -'+str(Vbias)+' V')
    return filteredSig

def plot_singlePulse(inputFile, Header, Vbias, cutoffFreq = 10e6, N = 1):    
    """
    A method to plot individual pulses from a CAEN wavedump file. 

    Parameters:
    ----------- 
    inputFile: string, input file-name
    N: int, number of events to be read 
    Header: 1 if input file has header, 0 otherwise 
    Vbias: biasing voltage 
    cutoffFreq: cutoff frequency in Hertz, 10 MHz is best for sipm paper
    
    Example use:
    -----------
    from signal_processing import plot_singlePulse
    vbias = 69; sipmN = '4'; source = 'dark'; cutoff = 10e6; 
    peak = 'single'; connection = 'series'; 
    if source == 'dark': 
        header = 1;
    else: 
        header = 0; 
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    plot_singlePulse(pathToData, header, vbias, cutoffFreq = cutoff);  
    """
    dataFile = wavedumpReader.DataFile(inputFile)
    dataFile.file.seek(0)
    plt.rcParams.update({'font.size': 12})
    plt.figure() 
    # Dark noise: 
    if Header == 1: 
        for i in range(N):
            header = fromfile(dataFile.file, dtype='I', count=6)
            if len(header) != 6:
                break
            eventSize = (header[0] - 24) // 2
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            fs = 256e6
            filteredSig = 1000*(filters.butter_lowpass_filter(trace,cutoffFreq,fs))
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.figure(); plt.plot(timeAx, filteredSig);
            # Unfiltered signal plot:
            #plt.plot(timeAx, trace); 
    # LED trigger:  
    else: 
        for i in range(N):
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=1024) 
            trace = trace - np.mean(trace[0:100]);
            trace = trace*2/(2**12);
            fs = 256e6
            filteredSig = 1000*(filters.butter_lowpass_filter(trace,cutoffFreq,fs))
            timeAx = linspace(0,4*len(trace),len(trace))
            # Filtered signal:
            plt.figure(); plt.plot(timeAx, filteredSig);
            # Unfiltered signal plot:
            #plt.plot(timeAx, trace);   
    del dataFile;
    plt.xlabel('Time (ns)'); plt.ylabel('Voltage (mV)')
    plt.title(str(cutoffFreq/(10**6))+' MHz BW filter, $V_{bias}$ = -'+str(Vbias)+' V')
    return filteredSig 
    
