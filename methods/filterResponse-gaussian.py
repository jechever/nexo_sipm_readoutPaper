from scipy import signal, ndimage
import matplotlib.pyplot as plt
import numpy as np

# Set peaking time and sigma input for Gaussian filter 
#peakingTime = 2000; sigma = peakingTime/15.; 

def get_peakingTime(peakingTime): 
    """
    A method to measure the peaking time response of a Gaussian filter. 
    Parameters: 
    ------------
    - peakingTime: Approximation for the desired peaking time. 

    Returns: 
    ------------
    pulse: average of 1000 pulses, rise time shown in plot

    Example use: 
    ------------
    from sipm_signalProcessing import getAvgPulse_gauss
    peMin = 0.021; peMax = 0.06; t_peak = 80
    vbias = 71; sipmN = '4'; source = 'dark'; connection = 'series'; 
    print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
    pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
    avg_pulse, t_rise = getAvgPulse_gauss(pathToData, peMin, peMax, sipmN, vbias, t_peak)

    Notes: 
    ------------
    This method only takes in data files ***with header***
    """
    sigma = peakingTime/15.; 
    # Generate delta function 
    nSamples, index = 1024, int(1024/2)
    imp = signal.unit_impulse(nSamples, index)
    time_ax_gauss = np.linspace(0,4*nSamples,nSamples)
    # Filter response: Gaussian
    response_gauss = ndimage.gaussian_filter1d(imp, sigma)
    plt.figure(); plt.plot(time_ax_gauss, imp, label = 'Dirac delta')
    #plt.plot(np.arange(-5000, 5000), response_gauss, label = 'Gauss filter response -Johny')
    #xAxis = (np.arange(-5000, 5000)+time_ax_gauss[int(detect_peaks.detect_peaks(adjusted_gauss, mph = 0.0002, mpd = 20000))])
    plt.plot(time_ax_gauss, response_gauss, label = 'Gaussian filter response -Johny')
    #plt.plot(time_ax_gauss, adjusted_gauss, label = 'Gauss filter response -Mietek')
    plt.margins(0.1, 0.1); plt.xlabel('Time [ns]'); plt.ylabel('Amplitude')
    plt.grid(True); plt.show(); plt.legend()
    #plt.xlim(-4000, 6000); 
    plt.ylim(-0.0001, max(response_gauss)*1.3);
    #plt.yscale('log');            
    for j in range(len(response_gauss)): 
        if response_gauss[j] >= 0.01*max(response_gauss):
            plt.axvline(x=time_ax_gauss[j], color='r',linestyle='--')
            minT = time_ax_gauss[j]
            break
    for j in range(len(response_gauss)):
        if response_gauss[j] == max(response_gauss):
            plt.axvline(x=time_ax_gauss[j], color='r',linestyle='--')
            maxT = time_ax_gauss[j]
            break            
    peakingTime =   maxT - minT;
    print('Peaking time: '+str(peakingTime))
    return peakingTime