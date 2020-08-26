"""
Single plot for all measurements of ENC vs T_peak with errorbars. Only 
measurements are shown, not the model. 
Author: J. Echevers, for the nEXO collaboration. For use outside of the 
collaboration contact author. 
"""

import matplotlib.pyplot as plt
import pandas as pd 

# Local path: 
# cd Documents/nEXO/sipm_analysis/code/scripts

## ################## ##
##     MAIN BODY      ##
## ################## ##

if __name__ == "__main__":
    
    # Read data file for measurements: 
    df_series_triang_6s = pd.read_csv('../data-output/enc-tp/triangular/enc-triangFilter-6s.csv')
    df_series_triang_4s = pd.read_csv('../data-output/enc-tp/triangular/enc-triangFilter-4s.csv')
    df_series_triang_2s = pd.read_csv('../data-output/enc-tp/triangular/enc-triangFilter-2s.csv')
    
    tpeak_series_triang = df_series_triang_6s['peaking time (ns)'].tolist()
    enc_series_triang_6sipm = df_series_triang_6s['ENC(s.p.e.)'].tolist()
    enc_series_triang_4sipm = df_series_triang_4s['ENC(s.p.e.)'].tolist()
    enc_series_triang_2sipm = df_series_triang_2s['ENC(s.p.e.)'].tolist()
    errors6 = df_series_triang_6s['error(s.p.e.)'].tolist()
    errors4 = df_series_triang_4s['error(s.p.e.)'].tolist()
    errors2 = df_series_triang_2s['error(s.p.e.)'].tolist()
    
    # Make ENC vs Tpeak plot: 
    plt.ion()
    plt.figure(); 
    plt.errorbar(tpeak_series_triang, enc_series_triang_6sipm, yerr=errors6, fmt='or', capsize = 2, capthick=2, label = '2s3p')
    plt.errorbar(tpeak_series_triang, enc_series_triang_4sipm, yerr=errors4, fmt='og', capsize = 2, capthick=2, label = '2s2p')
    plt.errorbar(tpeak_series_triang, enc_series_triang_2sipm, yerr=errors2, fmt='ok', capsize = 2, capthick=2, label = '2s')
    plt.title('Triangular filter, measurement, $V_{over}$ = 4.5 V') ; 
    plt.xlabel('$t_{peak}$ (s)'); plt.ylabel('enc ( s.p.e.)');
    plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()
    plt.savefig('../figures/enc-tp/triangular/enc-all-measurement.png')
    
    # Read model file: 
    df_series_triang = pd.read_csv('../noise-model/enc-triangFilter-noiseModel.csv')
    
    # Make lists
    tpeak_series_triang = df_series_triang['peaking time (ns)'].tolist()
    enc_series_triang_6sipm = df_series_triang['ENC-3p2s'].tolist()
    enc_series_triang_4sipm = df_series_triang['ENC-2p2s'].tolist()
    enc_series_triang_2sipm = df_series_triang['ENC-2s'].tolist()
    
    # Make ENC vs Tpeak plot: 
    plt.figure(); 
    plt.plot(tpeak_series_triang, enc_series_triang_6sipm, 'r', label = '2s3p')
    plt.plot(tpeak_series_triang, enc_series_triang_4sipm, 'g', label = '2s2p')
    plt.plot(tpeak_series_triang, enc_series_triang_2sipm, 'k', label = '2s')
    plt.title('Triangular filter, model, $V_{over}$ = 4.5 V') ; 
    plt.xlabel('$t_{peak}$ (s)'); plt.ylabel('enc ( s.p.e.)');
    plt.grid(); plt.show(); plt.yscale('log'); plt.xscale('log'); plt.legend()
    plt.savefig('../figures/enc-tp/triangular/enc-all-model.png')
    
    input('Press [ENTER] to finish.')
    
    