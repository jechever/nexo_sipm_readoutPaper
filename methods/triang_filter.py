# by M. Dabrowski 08/12/2020
import matplotlib.pyplot as plt
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


## ################## ##
##     MAIN BODY      ##
## ################## ##

# Uncomment for example: 

#"""
if __name__ == "__main__":

    N = 1024   # number of samples
    Tstep = 4e-9  # time-step

    # Generation of a dummy SiPM signal
    xsipm, sipm  = gen_sipm_sig(N, Tstep, 1e-7)

    xtf, tf0, tri0 = trig_filter(sipm, N, Tstep, 0.1e-6)
    xtf, tf1, tri1 = trig_filter(sipm, N, Tstep, 0.5e-6)
    xtf, tf2, tri2 = trig_filter(sipm, N, Tstep, 1e-6)

    plt.ion()
    
    plt.figure()
    plt.title('SiPM signal and filter weighting functions')
    plt.plot(xsipm, sipm, '-o', markersize=1.5, label='Dummy SiPM signal')
    plt.plot(xtf, tri0, '-o', markersize=1.5, label='Triang Filter - t$_p$=100 ns')
    plt.plot(xtf, tri1, '-o', markersize=1.5, label='Triang Filter - t$_p$=500 ns')
    plt.plot(xtf, tri2, '-o', markersize=1.5, label='Triang Filter - t$_p$=1 us')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel('time (s)'); plt.ylabel('a.u.'); 
    plt.ylim(0, 1.5);
    plt.draw()
    
    plt.figure()
    plt.title('Signal after triangular filtering')
    plt.plot(xtf, tf0, '-o', markersize=1.5, label='Filtered signal - t$_p$=100 ns')
    plt.plot(xtf, tf1, '-o', markersize=1.5, label='Filtered signal - t$_p$=500 ns')
    plt.plot(xtf, tf2, '-o', markersize=1.5, label='Filtered signal - t$_p$=1 us')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel('time (s)'); plt.ylabel('a.u.'); 
    plt.ylim(0, 35);
    plt.draw()
    
    #input("Press [enter] to FINISH.")
    #exit(0)
#"""