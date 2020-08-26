# by M. Dabrowski 08/12/2020

import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
from math import floor, exp

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
       
if __name__ == "__main__":
    
    df_pulse = pd.read_csv('sample-pulses/pulse_0.csv')
    pulse = df_pulse['pulse (mV)'].tolist()
    pulse = [-x for x in pulse]
    N = len(pulse)   # number of samples
    timeAx = np.linspace(0,4*len(pulse),len(pulse))
    Tstep = 4e-9  # time-step
    
    xtf, tf0, tri0 = trig_filter(pulse, N, Tstep, 0.1e-6)
    xtf, tf1, tri0 = trig_filter(pulse, N, Tstep, 0.2e-6)
    xtf, tf2, tri1 = trig_filter(pulse, N, Tstep, 0.5e-6)
    xtf, tf3, tri2 = trig_filter(pulse, N, Tstep, 1e-6)
    
    plt.ion()
    plt.figure()
    plt.title('Signal after triangular filtering')
    plt.plot(timeAx, tf0, '-o', markersize=1.5, label='Filtered signal - t$_p$=100 ns')
    plt.plot(timeAx, tf1, '-o', markersize=1.5, label='Filtered signal - t$_p$=200 ns')
    plt.plot(timeAx, tf2, '-o', markersize=1.5, label='Filtered signal - t$_p$=500 ns')
    plt.plot(timeAx, tf3, '-o', markersize=1.5, label='Filtered signal - t$_p$=1 us')
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel('time (ns)'); plt.ylabel('a.u.'); 
    plt.ylim(-0.2, 2);
    plt.draw()
    
    plt.figure()
    plt.title('Unfiltered signal')
    plt.plot(xtf, pulse, '-o', markersize=1.5)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel('time (s)'); plt.ylabel('a.u.'); 
    #plt.ylim(0, 35);
    plt.draw()
