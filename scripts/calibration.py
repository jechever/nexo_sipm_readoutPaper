"""
Quick check script for CAEN Wavedump waveforms
"""
import sys

# Set path to SiPM methods
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/Users/johny/Documents/nEXO/sipm_analysis/code/methods')
from signal_processing import get_amplitudes
from matplotlib.ticker import FormatStrFormatter
from plotting import plot_rawPulses
import matplotlib.pyplot as plt
 
# Local path: 
# cd Documents/nEXO/sipm_analysis/code

# Set global variables: 
Npulses = 20000;       # analyze the first 20000 pulses from the data set
vbias = 68;            # biasing voltage
connection = 'series'; #'parallel' or 'series'
sipmN = '6';           # Number of SiPMs
source = 'dark';       # 'dark' or 'led', depending on how data was taken
# Set path to data file: 
pathToData = '../../sipm_data/calibration/data_10PF_0.01Vpp_amp.dat'
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')

if __name__ == "__main__":
    
    plt.ion()
    
    # Get amplitudes to generate histogram 
    #amplitudes, integrals = get_amplitudes(pathToData, Npulses, header, vbias, fullWindow = True); 
    amplitudes, integrals = get_amplitudes(pathToData, Npulses, 1, vbias, t_peak = 200, filtering = 'triangular', plotting = 1);
    
    #plot_rawPulses(pathToData, header, vbias) 
    plot_rawPulses(pathToData, 1, vbias) 
    
    # Plot PE spectrum:
    #lower, upper = -5, 200; 
    nbins = 100; 
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.rcParams.update({'font.size': 12})
    y,x,_= plt.hist(amplitudes, nbins, histtype = 'step', density=True, label = 'V$_{bias}$ = -'+str(vbias)+' V');
    plt.xlabel('Pulse Amplitude (V)'); plt.ylabel('Events/bin')
    plt.legend(); plt.show()

    input('Press [ENTER] to finish.')