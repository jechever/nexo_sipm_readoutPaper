"""
Quick check script for CAEN Wavedump waveforms
"""
import sys

# Set path to SiPM methods
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../methods')
from signal_processing import get_amplitudes
from matplotlib.ticker import FormatStrFormatter
from plotting import plot_rawPulses
import matplotlib.pyplot as plt
 
# Local path: 
# cd Documents/nEXO/sipm_analysis/code/scripts

# Set global variables: 
Npulses = 20000;       # analyze the first 20000 pulses from the data set
peakingTime = 100; 

# Set path to data file: 
pathToData = '../../sipm_data/calibration/data_10PF_0.01Vpp_amp.dat'

if __name__ == "__main__":
    
    plt.ion()
    
    # Get amplitudes to generate histogram  
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, 1, t_peak = peakingTime, plotting = 1);
    
    #plot_rawPulses(pathToData, header, vbias) 
    plot_rawPulses(pathToData, 1) 
    
    # Plot PE spectrum:
    #lower, upper = -5, 200; 
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', density=True);
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.show()

    input('Press [ENTER] to finish.')