import numpy as np
from sipm_signalProcessing import getAvgPulse
Tlist = []

peMin = 0.021; peMax = 0.06; 
vbias = 71; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias)
Tlist.append(Trise)

peMin = 0.018; peMax = 0.054; 
vbias = 70; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias)
Tlist.append(Trise)

peMin = 0.018; peMax = 0.051; 
vbias = 69; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias)
Tlist.append(Trise)

peMin = 0.018; peMax = 0.046; 
vbias = 68; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias)
Tlist.append(Trise)

peMin = 0.015; peMax = 0.041; 
vbias = 67; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias)
Tlist.append(Trise)

peMin = 0.012; peMax = 0.037; 
vbias = 66; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias)
Tlist.append(Trise)

peMin = 0.011; peMax = 0.03; 
vbias = 65; sipmN = '4'; source = 'dark'; connection = 'series'; 
print('Analysing '+str(sipmN)+'-sipms biased at '+str(vbias)+'V')
pathToData = '../sipm_data/'+sipmN+'sipm_testing/'+connection+'_test/wave0_'+str(vbias)+'V_'+source+'.dat';   
pulse, Trise = getAvgPulse(pathToData, peMin, peMax, sipmN, vbias)
Tlist.append(Trise)

trise, sigma = np.mean(Tlist), np.std(Tlist)
print('rise time is: '+str(trise)+'+/-'+str(sigma)+' ns')