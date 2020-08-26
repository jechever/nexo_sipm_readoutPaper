from numpy import fromfile, dtype, linspace, mean
from scipy import ndimage
#from numba import jit
import matplotlib.pyplot as plt
import numpy as np 
import wavedumpReader
import detect_peaks
import filters     
    
#-------------------------------------------------------------------------------
#--------------------------- Analysis methods ----------------------------------
#-------------------------------------------------------------------------------

class DataFile:
    def __init__(self, fileName):
        """
        Initializes the dataFile instance to include the fileName, access time,
        and the number of boards in the file. Also opens the file for reading.
        The file must have been collected with the OUTPUT_FILE_HEADER option
        set to YES. Otherwise, the reader will fail.
        """
        self.fileName = path.abspath(fileName)
        self.file = open(self.fileName, 'rb')
        self.recordLen = 0
        self.oldTimeTag = 0.
        self.timeTagRollover = 0
        self.boardId = 0
        
        self.T = T        
        self.fileName = fileName
        self.PEmin = PEmin, 
        self.PEmax = PEmax
        self.Nsipms = , 
        
        Vbias, 
        min_rise = 0.01, 
        cutoffFreq = 10e6, 
        connection = 'series', 
        Header = 0, 
        plotPulses = 0
        
    
def getAvgPulse_gauss(fileName, PEmin, PEmax, Nsipms, Vbias, t_peak, min_rise = 0.01, connection = 'series', Header = 0, plotPulses = 0):
    """
    A method to extract the pulse rise time (Tau) in units of ns from 
    CAEN wavedump files. 

    Parameters: 
    ------------
    - fileName: input file 
    - PEmin: lower limit for SPE 
    - PEmax: lower limit for SPE
    - Nsipms: number of SiPMs
    - Vbias: biasing voltage

    Returns: 
    ------------
    pulse: average of 1000 pulses, rise time shown in plot

    Example use: 
    ------------
    from signal_processing import getAvgPulse_gauss
    peMin = 0.021; peMax = 0.06; t_peak = 80
    vbias = 71; sipmN = '4'; source = 'dark'; connection = 'series'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
    avg_pulse, t_rise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, t_peak)

    Notes: 
    ------------
    This method only takes in data files ***with header***
    """
    dataFile = wavedumpReader.DataFile(fileName)
    dataFile.file.seek(0) 
    i = 0; pulses = [];  
    sigma = t_peak/15
    if connection == 'series':
        Vover = (Vbias - 59)/2
    if Header == 1: 
        while True:#
            header = fromfile(dataFile.file, dtype='I', count=6)
            if len(header) != 6:
                break
            eventSize = (header[0] - 24) // 2
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize)
            trace = trace - mean(trace[0:100]);
            trace = trace*2/(2**12); 
            
            T = linspace(0,4*len(trace),len(trace))
            filteredSig = ndimage.gaussian_filter1d(trace, sigma)#
            ind_valleys = detect_peaks.detect_peaks(filteredSig, mph = PEmin, valley=True)
            #---------------------
            if abs(min(filteredSig)) > abs(PEmin) and abs(PEmax) > abs(min(filteredSig)):
                if len(ind_valleys) == 1:
                    pulses.append(filteredSig);
                    i += 1;
                    if i == 1000: 
                        break
        pulse = 1000*(np.mean(pulses, axis=0))
        if plotPulses != 0:
            plt.figure(); timeAx = T;             
            for j in range(len(pulse)): 
                if pulse[j] <= min_rise*min(pulse):
                    plt.axvline(x=j*4, color='r',linestyle='--')
                    minT = j*4
                    break
            for j in range(len(pulse)):
                if pulse[j] == min(pulse):
                    plt.axvline(x=j*4, color='r',linestyle='--')
                    maxT = j*4
                    break
        elif plotPulses == 0:     
            for j in range(len(pulse)): 
                if pulse[j] <= min_rise*min(pulse):
                    minT = j*4
                    break
            for j in range(len(pulse)):
                if pulse[j] == min(pulse):
                    maxT = j*4
                    break
    elif Header == 0: 
        while True:
            eventSize = 1024
            trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize)
            trace = trace - mean(trace[0:100]);
            trace = trace*2/(2**12);
            T = linspace(0,4*len(trace),len(trace))
            filteredSig = ndimage.gaussian_filter1d(trace, sigma)#
            ind_valleys = detect_peaks.detect_peaks(filteredSig, mph = PEmin, valley=True)
            #---------------------
            if abs(min(filteredSig)) > abs(PEmin) and abs(PEmax) > abs(min(filteredSig)):
                if len(ind_valleys) == 1:
                    pulses.append(filteredSig);
                    i += 1;
                    if i == 1000: 
                        break
        pulse = 1000*(np.mean(pulses, axis=0))
        if plotPulses != 0:
            plt.figure(); timeAx = T; 
            for j in range(len(pulse)): 
                if pulse[j] <= min_rise*min(pulse):
                    plt.axvline(x=j*4, color='r',linestyle='--')
                    minT = j*4
                    break
            for j in range(len(pulse)):
                if pulse[j] == min(pulse):
                    plt.axvline(x=j*4, color='r',linestyle='--')
                    maxT = j*4
                    break
        elif plotPulses == 0:     
            for j in range(len(pulse)): 
                if pulse[j] <= min_rise*min(pulse):
                    minT = j*4
                    break
            for j in range(len(pulse)):
                if pulse[j] == min(pulse):
                    maxT = j*4
                    break
    tau = maxT-minT
    if plotPulses != 0:
        plt.plot(timeAx, pulse, label='$T_{peak} = $'+str(tau)+' (ns)'); 
        plt.title('$V_{over} = $'+str(Vover)+' V, Gauss filter')
        plt.xlabel('time (ns)'); plt.ylabel('Voltage (mV)'); 
        plt.show(); plt.legend(); 
    print('Peaking time: ', tau,' ns')
    return pulse, tau


def get_amplitudes(inputFile, N, Header, Vbias, t_peak, cutoffFreq = 10e6, fullWindow = True, filtering = 'triangular', plotting = 0):    
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
    from  import get_amplitudes, plot_peSpectrum 
    Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series';  
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias);
    nbins = 200; lower = 0.; upper = 80; 
    plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
    """
    dataFile = wavedumpReader.DataFile(inputFile)
    dataFile.file.seek(0)
    plt.rcParams.update({'font.size': 12})
    amplitude_list = []; integral_list = []  
    t_peak *= 1e-9; 
    #minR = 470; maxR = 515;
    if plotting != 0: 
            plt.figure();
    #for i in range(N):
    for i in range(N):
        if Header == 1: 
            header = fromfile(dataFile.file, dtype='I', count=6)
            if len(header) != 6:
                break
            eventSize = (header[0] - 24) // 2
            minR = 470; maxR = 515;
        else: 
            eventSize = 1024
            minR = 500; maxR = 610;
        minR = 470; maxR = 515;
        trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
        trace = trace - np.mean(trace[0:100]);
        trace = trace*2/(2**12); trace = [-x for x in trace]
        if filtering == 'bw': 
            fs = 256e6
            pulse = 1000*(filters.butter_lowpass_filter(trace,cutoffFreq,fs))
        elif filtering == 'triangular': 
            N = eventSize; Tstep = 4e-9; 
            xtf, pulse, tri0 = filters.trig_filter(trace, N, Tstep, t_peak)
        else: 
            sigma = t_peak/15.
            pulse = 1000*(ndimage.gaussian_filter1d(trace, sigma))#
        timeAx = linspace(0,4*len(trace),len(trace))
        if fullWindow == True: 
            amplitude = max(pulse); amplitude_list.append(amplitude)
            integral = np.trapz(pulse, x=timeAx); integral_list.append(integral)
        else: 
            amplitude = max(pulse[minR:maxR]); amplitude_list.append(amplitude)
            integral = np.trapz(pulse[minR:maxR], x=timeAx[minR:maxR]); integral_list.append(integral)
        if plotting != 0 and i < 100:                         
            #plt.plot(timeAx, trace)
            plt.plot(timeAx, pulse)
            plt.xlabel('Time (ns)'); plt.ylabel('Voltage (mV)')
    del dataFile;
    return amplitude_list, integral_list

def get_amplitudes_noAP(inputFile, N, Header, Vbias, peMin, cutoffFreq = 10e6, t_peak = 100, filtering = 'gauss', fullWindow = True, plotting = 0):    
    """
    A method to obtain a list of amplitudes and integrals from FBK SiPMs, 
    formatted for CAEN Wavedump readout. 

    Parameters:
    ----------- 
    inputFile: string, input file-name
    N: int, number of events to be read 
    Header: 1 if input file has header, 0 otherwise 
    Vbias: biasing voltage 
    cutoffFreq: cutoff frequency in Hertz, 10 MHz is best for sipm paper

    Example use:
    -----------
    from signal_processing import get_amplitudes, plot_peSpectrum 
    Npulses = 20000; vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
    if source == 'dark': 
        minR = 470; maxR = 515;
        header = 1;
    else: 
        minR = 560; maxR = 610;
        header = 0; 
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    amplitudes, integrals = get_amplitudesFullWindow(pathToData, Npulses, header, vbias);
    nbins = 200; lower = 0.; upper = 80; 
    plot_peSpectrum(amplitudes, nbins, lower, upper, vbias, connection, sipmN);
    """
    dataFile = wavedumpReader.DataFile(inputFile)
    dataFile.file.seek(0)
    plt.rcParams.update({'font.size': 12})
    amplitude_list = []; integral_list = []; 
    if filtering == 'bw': 
        # Dark noise: 
        if Header == 1: 
            minR = 470; maxR = 515; i = 0; 
            while True: 
                header = fromfile(dataFile.file, dtype='I', count=6)
                if len(header) != 6:
                    break
                eventSize = (header[0] - 24) // 2
                trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
                trace = trace - np.mean(trace[0:100]);
                trace = trace*2/(2**12);
                fs = 256e6
                filteredSig = 1000*(filters.butter_lowpass_filter(trace,cutoffFreq,fs))
                pulse = [-x for x in filteredSig]
                timeAx = linspace(0,4*len(trace),len(trace))
                #if i == 1:
                    #plt.figure(); plt.plot(timeAx, pulse)
                ind_valleys = detect_peaks.detect_peaks(pulse, mph = peMin*1000)
                if len(ind_valleys) == 1 and i < 20000:
                    if fullWindow == True: 
                        amplitude = max(pulse); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse, x=timeAx); integral_list.append(integral)
                    elif fullWindow == False: 
                        amplitude = max(pulse[minR:maxR]); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse[minR:maxR], x=timeAx[minR:maxR]); integral_list.append(integral)
                    i += 1
        # LED trigger:  
        elif Header == 0: 
            minR = 560; maxR = 610;
            while True: 
                trace = fromfile(dataFile.file, dtype=dtype('<H'), count=1024) 
                trace = trace - np.mean(trace[0:100]);
                trace = trace*2/(2**12);
                fs = 256e6
                filteredSig = 1000*(filters.butter_lowpass_filter(trace,cutoffFreq,fs))
                pulse = [-x for x in filteredSig]
                timeAx = linspace(0,4*len(trace),len(trace))
                #if i == 1:
                    #plt.figure(); plt.plot(timeAx, pulse)
                ind_valleys = detect_peaks.detect_peaks(pulse, mph = peMin*1000)
                if len(ind_valleys) == 1 and i < 20000:
                    if fullWindow == True: 
                        amplitude = max(pulse); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse, x=timeAx); integral_list.append(integral)
                    elif fullWindow == False: 
                        amplitude = max(pulse[minR:maxR]); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse[minR:maxR], x=timeAx[minR:maxR]); integral_list.append(integral)
                    i += 1
    elif filtering == 'gauss': 
        sigma = t_peak/15
        # Dark noise: 
        if Header == 1: 
            minR = 470; maxR = 515; i = 0; 
            if plotting != 0:
                plt.figure()
            while True: 
                header = fromfile(dataFile.file, dtype='I', count=6)
                if len(header) != 6:
                    break
                eventSize = (header[0] - 24) // 2
                trace = fromfile(dataFile.file, dtype=dtype('<H'), count=eventSize) 
                trace = trace - np.mean(trace[0:100]);
                trace = trace*2/(2**12);
                ind_valleys = detect_peaks.detect_peaks(trace, mph = peMin, valley=True)
                filteredSig = filteredSig = 1000*(ndimage.gaussian_filter1d(trace, sigma))
                pulse = [-x for x in filteredSig]
                timeAx = linspace(0,4*len(trace),len(trace))
                if len(ind_valleys) == 1 and i < 20000 and abs(min(trace[0:445])) < peMin:
                    # Sanity check:
                    if plotting != 0 and i < 100:                         
                        trace = [-x*1000 for x in trace]
                        #plt.plot(timeAx, trace)
                        plt.plot(timeAx, pulse)
                        plt.xlabel('Time (ns)'); plt.ylabel('Voltage (mV)')
                    if fullWindow == True: 
                        amplitude = max(pulse); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse, x=timeAx); integral_list.append(integral)
                    elif fullWindow == False: 
                        amplitude = max(pulse[minR:maxR]); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse[minR:maxR], x=timeAx[minR:maxR]); integral_list.append(integral)
                    i += 1
                if i == 2000: 
                    break 
        # LED trigger:  
        elif Header == 0: 
            minR = 525; maxR = 610; i = 0;
            if plotting != 0:
                plt.figure()
            while True: 
                trace = fromfile(dataFile.file, dtype=dtype('<H'), count=1024) 
                trace = trace - np.mean(trace[0:100]);
                trace = trace*2/(2**12);
                ind_valleys = detect_peaks.detect_peaks(trace, mph = peMin, valley=True)
                filteredSig = filteredSig = 1000*(ndimage.gaussian_filter1d(trace, sigma))
                pulse = [-x for x in filteredSig]
                timeAx = linspace(0,4*len(trace),len(trace))
                if len(ind_valleys) == 1 and i < 2000 and abs(min(trace[0:445])) < peMin:
                    # Sanity check:
                    if plotting != 0 and i < 100: 
                        trace = [-x*1000 for x in trace]
                        #plt.plot(timeAx, trace)
                        plt.plot(timeAx, pulse)
                        plt.xlabel('Time (ns)'); plt.ylabel('Voltage (mV)')
                    if fullWindow == True: 
                        amplitude = max(pulse); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse, x=timeAx); integral_list.append(integral)
                    elif fullWindow == False: 
                        amplitude = max(pulse[minR:maxR]); amplitude_list.append(amplitude)
                        integral = np.trapz(pulse[minR:maxR], x=timeAx[minR:maxR]); integral_list.append(integral)
                    i += 1
                if i == 2000: 
                    break 
    del dataFile;
    return amplitude_list, integral_list
    
