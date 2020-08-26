"""
Methods for implementing digital filters in SiPM 2020 paper  
"""

from scipy.signal import butter, lfilter, filtfilt
from math import floor, exp
import numpy as np

#############################
##### Triangular filter #####
#############################

def gen_triang_wf(NS, TS, tp):
    xtri = [ (x+1)*TS for x in range(-floor(NS/2), floor(NS/2))]    
    tri_ = []
    tri_.append(1)    
    for i in range(1,floor(NS/2)):
        if i*TS<tp:
            y = (tp - i*TS)/tp
        else:
            y = 0
        tri_.append(y)
    tri = []
    for y_ in np.flip(tri_[1:],0):
        tri.append(y_)
    for y_ in tri_:
        tri.append(y_)
    tri.append(0)
    #tri = tri / np.sum(tri)    
    return xtri, tri

# Dummy time-domain SiPM signal
def gen_sipm_sig(NS, TS, tau):
    X = [ (x+1)*TS for x in range(-floor(NS/2), floor(NS/2))]
    sig = []   
    for x in X:
        if x < 0:
            sig.append(0)
        else:
            sig.append(exp(-x/tau))
    return X, sig

# Returns X and Y values
def trig_filter(input, NS, TS, tp):
    # NS - number of samples
    # TS - time step between samples
    # tp - filter peaking time    
    xtri, tri = gen_triang_wf(NS, TS, tp)
    return xtri, np.convolve(tri, input, mode='same'), tri

#############################
#### Butterworth filter #####
#############################

def butter_lowpass_filter(data, cutoff, fs, order=5):
    """
    Parameters:
    ----------- 
        data: signal, can be a numpy array or a list
        cutoff: cutoff frequency (Hertz) 
        fs: sampling frequency (samples/s)
        order: BW filter order to use
        
    Returns: 
    -----------
        y: filtered signal 
        
    Example use: 
    -----------
        # Filter a signal with sampling freq (fs) of 256 Mega samples/s
        # using a 10 MHz lowpass BW filter. 
        import BWfilter
        cutoffFreq = 10e6; fs = 256e6
        filtered_signal = BWfilter.butter_lowpass_filter(trace,cutoffFreq,fs)
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
    
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def rms(array):
    from numpy import mean, sqrt, square
    rms = sqrt(mean(square(array)))
    return rms 

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
    
def butter_highpass(cutoff, fs, order=5): 
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

