"""
Script for getting IVs for 2019 large-area SiPM paper
"""

import wavedumpReader
import sipm_dataAnalysis
import detect_peaks
import matplotlib.pyplot as plt


# ------------------------------------------
# ---- Make IV curves for carriers 1&2 -----
# ------------------------------------------

Vbias = [33, 33.5, 34, 34.5, 35, 35.5, 36]

#Carrier 1 
i_s6 = [64E-9, 70E-9, 0.14E-6, 0.3E-6, 0.54E-6, 0.85E-6, 1.2E-6]
i_s5 = [74E-9, 78E-9, 0.14E-6, 0.29E-6, 0.52E-6, 0.83E-6, 1.2E-6]
i_s4 = [78E-9, 85E-9, 0.15E-6, 0.34E-6, 0.62E-6, 0.98E-6, 1.4E-6]
i_s1 = [89E-9, 96E-9, 0.20E-6, 0.43E-6, 0.8E-6, 1.2E-6, 1.8E-6]
i_s2 = [90E-9, 95E-9, 0.15E-6, 0.31E-6, 0.57E-6, 0.92E-6, 1.3E-6]
plt.figure(); 
plt.plot(Vbias, i_s1, label = 's1'); 
plt.plot(Vbias, i_s2, label = 's2'); 
plt.plot(Vbias, i_s4, label = 's4'); 
plt.plot(Vbias, i_s5, label = 's5'); 
plt.plot(Vbias, i_s6, label = 's6'); 
plt.xlabel('|V$_{bias}$| (V)'); plt.ylabel('|I| (A)')
plt.grid(); plt.show(); plt.yscale('log'); 
plt.legend()

#Carrier 2
i2_s1 = [90E-9, 97E-9, 0.17E-6, 0.34E-6, 0.6E-6, 0.96E-6, 1.4E-6]
i2_s2 = [91E-9, 98E-9, 0.17E-6, 0.35E-6, 0.62E-6, 0.99E-6, 1.47E-6]
i2_s4 = [97E-9, 99E-9, 0.16E-6, 0.33E-6, 0.6E-6, 0.97E-6, 1.4E-6]
i2_s5 = [92E-9, 98E-9, 0.16E-6, 0.33E-6, 0.57E-6, 0.9E-6, 1.23E-6]
plt.figure(); 
plt.title('carrier 2')
plt.plot(Vbias, i2_s1, label = 's1'); 
plt.plot(Vbias, i2_s2, label = 's2'); 
plt.plot(Vbias, i2_s4, label = 's4'); 
plt.plot(Vbias, i2_s5, label = 's5'); 
plt.xlabel('|V$_{bias}$| (V)'); plt.ylabel('|I| (A)')
plt.grid(); plt.show(); plt.yscale('log'); 
plt.legend()

#Both carriers
plt.figure(); 
plt.plot(Vbias, i2_s1, 'k:', label = 's1, carrier 2'); 
plt.plot(Vbias, i2_s2, 'k:', label = 's2, carrier 2'); 
plt.plot(Vbias, i2_s4, 'k:', label = 's4, carrier 2'); 
plt.plot(Vbias, i2_s5, 'k:', label = 's5, carrier 2'); 
plt.plot(Vbias, i_s1, 'b:', label = 's1, carrier 1'); 
plt.plot(Vbias, i_s2, 'b:', label = 's2, carrier 1'); 
plt.plot(Vbias, i_s4, 'b:', label = 's4, carrier 1'); 
plt.plot(Vbias, i_s5, 'b:', label = 's5, carrier 1');
plt.plot(Vbias, i_s6, 'b:',label = 's6, carrier 1'); 
plt.xlabel('|V$_{bias}$| (V)'); plt.ylabel('|I| (A)')
plt.grid(); plt.show(); plt.yscale('log'); 
plt.legend()
