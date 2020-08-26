"""
Script for Dark Noise analysis for 2020 SiPM paper 
"""
from nuissanceParams_methods import getDarkNoise_test, getTimeArray 

# SET PATH 
# cd Documents/nEXO/sipm_analysis/code

peMin = 0.018; peMax = 0.054; 
vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
time = getTimeArray(pathToData, peMin, peMax, sipmN, vbias)
getDarkNoise_test(time, sipmN, vbias)
#plt.yscale('log')
