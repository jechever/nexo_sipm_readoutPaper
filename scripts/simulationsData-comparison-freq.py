import matplotlib.pyplot as plt
import pandas as pd 


df_m = pd.read_csv('200605_enc_spe_in_freq_butt_2_in_series.csv')
df_johny = pd.read_csv('enc-tpeak-AP-4.5-OV.csv')

freq_mietek = df_m['freq'].tolist()
enc_mietek = df_m['enc_4sipm'].tolist()
enc_johny = df_johny['resolution (s.p.e.)'].tolist()
freq_johny = df_johny['frequency cut (MHz)'].tolist()

plt.figure(); 
plt.plot(freq_johny, enc_johny, label = 'Data ' ); 
plt.plot(freq_mietek, enc_mietek, label = 'simulations from Mietek' );
plt.yscale('log'); plt.xscale('log')
plt.grid(True); plt.show(); plt.legend()