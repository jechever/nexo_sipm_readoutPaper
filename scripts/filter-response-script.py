from scipy import signal
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import detect_peaks 
import BWfilter

peakingTime = 1000; sigma = peakingTime/15.; 

# Read data from Mietek: BW filter
df_bw_response = pd.read_csv('resp_butterworth_1Mhz.csv')
df_bw_response.time *=1e9
mietek = df_bw_response[' V'].tolist()
time_ax = df_bw_response['time'].tolist()

# Read data from Mietek: Gauss filter
df_gauss_response = pd.read_csv('simulations_mietek/5thorder_semigauss_shaper_1us.csv')
df_gauss_response.time *=1e9
pulse_gauss_mietek = df_gauss_response['V'].tolist()
time_ax_gauss = df_gauss_response['time'].tolist()

# Generate delta function 
#dirac = signal.unit_impulse(500015)
index = int(10000/2+190)
imp = signal.unit_impulse(10000, index)

# Filter response: BW
response_bw = BWfilter.butter_lowpass_filter(imp,0.99e6,1000e6)
alpha = max(response_bw)/max(mietek) 
adjusted_bw = [x* alpha for x in mietek]

# Filter response: Gaussian
response_gauss = gaussian_filter1d(imp, sigma)
alpha_gauss = max(response_gauss)/max(pulse_gauss_mietek) 
adjusted_gauss = [x*alpha_gauss for x in pulse_gauss_mietek]

# Plot signals for BW filter 
plt.figure()
plt.plot(np.arange(-5000, 5000), imp, label = 'Dirac delta')
plt.plot(np.arange(-5000, 5000), response_bw, label = 'BW filter response -Johny')
plt.plot(time_ax, adjusted_bw, label = 'BW filter response -Mietek')
plt.margins(0.1, 0.1); plt.xlabel('Time [ns]'); plt.ylabel('Amplitude')
plt.grid(True); plt.show(); plt.legend()
plt.xlim(0, 3000); plt.ylim(-0.002, 0.004);
#plt.yscale('log');

# Plot signals for Gaussian filter
plt.figure()
plt.plot(np.arange(-5000, 5000), imp, label = 'Dirac delta')
#plt.plot(np.arange(-5000, 5000), response_gauss, label = 'Gauss filter response -Johny')
xAxis = (np.arange(-5000, 5000)+time_ax_gauss[int(detect_peaks.detect_peaks(adjusted_gauss, mph = 0.0002, mpd = 20000))])
plt.plot(xAxis, response_gauss, label = 'Gauss filter response -Johny')
plt.plot(time_ax_gauss, adjusted_gauss, label = 'Gauss filter response -Mietek')
plt.margins(0.1, 0.1); plt.xlabel('Time [ns]'); plt.ylabel('Amplitude')
plt.grid(True); plt.show(); plt.legend()
plt.xlim(-4000, 6000); plt.ylim(-0.0001, max(response_gauss)*1.3);
#plt.yscale('log');
        
for j in range(len(response_gauss)): 
    if response_gauss[j] >= 0.01*max(response_gauss):
        plt.axvline(x=xAxis[j], color='r',linestyle='--')
        minT = xAxis[j]
        break
for j in range(len(response_gauss)):
    if response_gauss[j] == max(response_gauss):
        plt.axvline(x=xAxis[j], color='r',linestyle='--')
        maxT = xAxis[j]
        break
        
peakingTime =   maxT - minT;
print('Peaking time: '+str(peakingTime))