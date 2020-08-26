from pylab import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np 

#------------------------------------
#------------Fit Functions-----------
#------------------------------------

# Standard Gaussian function
def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)
    
# A two-peak gaussian fit   
def diModal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

# A three-peak gaussian     
def triModal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)

# A four-peak gaussian     
def fourModal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3,mu4,sigma4,A4):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)+gauss(x,mu4,sigma4,A4)

#------------------------------------
#----------Fitting Methods-----------
#------------------------------------
   
def singleFitter(xFit, yFit, expected, plotting = 1): 
    """
    A method that fits a single Gaussian peak to input data.
    
    Parameters:
    ----------- 
    x: list or array type, array from histo 
    y: list or array type, array from histo 
    expected: list or array type, expected fit values for mu, sigma, A
    
    Returns: 
    ------------
    params: parameters corresponding to fit, each row is as follows: (mu, sigma, A) 
    plots fit on top of distribution
    sigma: standard error of parameters
     
    Example use:
    ------------ 
    import sipm_dataAnalysis
    y_34_d4, x_34_d4 = sipm_dataAnalysis.phsPlotter(filtPHS_34V_d4, 100, -0.184,-0.026,34);
    delta = [];
    delt, sig = singleFitter(x_34_d4[88:100],y_34_d4[88:100],x_34_d4,y_34_d4,(-0.04,0.01,540)); delta.append(delt[0])
    sigma = delt[1];
    delt, sig = singleFitter(x_34_d4[63:76],y_34_d4[63:76],x_34_d4,y_34_d4,(-0.08,0.01,250)); delta.append(delt[0])
    delt, sig = singleFitter(x_34_d4[40:51],y_34_d4[40:51],x_34_d4,y_34_d4,(-0.12,0.01,125)); delta.append(delt[0])
    print 'resolution is: ', sigma/mean([abs(delta[2]-delta[1]),abs(delta[1]-delta[0])])
    print 'peak to peak is: ', mean([abs(delta[2]-delta[1]),abs(delta[1]-delta[0])])
    """
    if len(yFit)>len(xFit):
        yFit = yFit[0:len(xFit)];
    else:
        xFit = xFit[0:len(yFit)];
    params,cov=curve_fit(gauss,xFit,yFit,expected);
    sigma=np.sqrt(np.diag(cov));
    x = np.linspace(xFit[0],xFit[-1],100)
    plt.plot(x,gauss(x,*params),color='red',lw=0.8); plt.legend();
    if plotting == 0: 
        plt.close()
    return params, sigma 

def dualGaussFit(xFit, yFit, xFull, yFull, expected):
    """
    A method to fit a function that is a sum of two Gaussians 
    
    Parameters: 
    -----------
    xFit: numpy array containing x-values of PHS, obtained from phsPlotter method
    yFit: numpy array containing y-values of PHS, obtained from phsPlotter method
    expected: numpy array containing rough fit values for two-Gaussian sum 
    
    Returns: 
    -----------
    params: parameters from fits for dual Gaussian function, each row is as follows: (mu, sigma, A) 
    params2: parameters from fits for Gaussian function corresponding to
             aferpulsing shoulder, each row is as follows: (mu, sigma, A) 
    plots fit on top of distribution, red dashed dual Guassian, green Gauss is afterpulsing shoulder 
    
    Example use: 
    -----------
    filtPHS34V_s6, singlePulse34V_s6, afterPulse35V_s6 = filterSpectraBin('jan_testing/ledOFF/withHeader/wave0-34V-sample6.dat', 20000, 0.03, 0.029, 470, 560, 1);
    y34_s6, x34_s6 = phsPlotter(singlePulse34V_s6, 100, -0.22, -0.029, 34); 
    dualGaussFit(x34_s6[73:99],y34_s6[73:99],x34_s6, y34_s6,[-0.05,0.01,1200,-0.06,0.01,310]);
    dualGaussFit(x34_s6[50:73],y34_s6[50:73],x34_s6, y34_s6,[-0.09,0.01,220,-0.11,0.01,130]);
    """ 
    if len(yFit)>len(xFit):
        yFit = yFit[0:len(xFit)];
    else:
        xFit = xFit[0:len(yFit)];
    params,cov=curve_fit(diModal,xFit,yFit,expected);
    params1 = [params[0], params[1], params[2]]
    params2 = [params[3], params[4], params[5]]
    sigma = params[1]; mu = params[0]
    plt.plot(xFull,gauss(xFull,*params1),color='green',lw=0.8, linestyle='dashed'); plt.legend();#No label
    plt.plot(xFull,gauss(xFull,*params2),color='green',lw=0.8, linestyle='dashed'); plt.legend();#No label
    plt.plot(xFull,diModal(xFull,*params),color='red',lw=0.9); plt.legend();#No label
    print('Sigma, and mu for SPE are: ')
    print(sigma, mu)
    #print('Parameters for afterpulsing Gaussian are: ')
    #print(params2)
    return params

# Four peak fitter
def fourFitter(x,y,expected): 
    if len(y)>len(x):
        y = y[0:len(x)];
    else:
        x = x[0:len(y)];
    params,cov=curve_fit(fourModal,x,y,expected);
    sigma=np.sqrt(diag(cov));
    plt.plot(x,fourModal(x,*params),color='red',lw=2,label='fit'); plt.legend();
    print(params);#,'\n',sigma)
    return params#, sigma 

# Four peak fitter
def triFitter(x,y,expected): 
    import numpy as np
    if len(y)>len(x):
        y = y[0:len(x)];
    else:
        x = x[0:len(y)];
    params,cov=curve_fit(triModal,x,y,expected);
    sigma=np.sqrt(diag(cov));
    plt.plot(x,triModal(x,*params),color='red',lw=2,label='fit'); plt.legend();
    print(params);#,'\n',sigma)
    return params#, sigma 
 

def resolution(params):
    """
    Numbers refer to parameters' array (mu, sigma, A,...):
        0 1 2  #1st photopeak
        3 4 5  #2nd...
        6 7 8 
        9 10 11
    """
    delta = (abs(params[0]-params[3])+abs(params[3]-params[6])+abs(params[6]-params[9]))/3; #FOR 4-fitter
    #delta = (abs(params[3]-params[6])+abs(params[6]-params[9]))/2; #FOR 3-fitter
    res = params[1]/delta; #For plots with NO pedestal, 1st peak -> 1st photopeak
    #res = params[4]/delta; #For plots with pedestal
    print(res);
    print(delta);
    return res, delta

def resolution_3peak(params):
    delta = (abs(params[0]-params[3])+abs(params[3]-params[6]))/3;
    res = params[4]/delta;
    print(res);
    print(delta);
    return res, delta

#plt.rcParams.update({'font.size': 14}) #To adjust plot fonts
"""     
**Working fit (using fourModal): 

#FIRST do (after changing the input file name in sipm_dataAnalysis):
     
%run "/Users/johny/Documents/nEXO/sipm_dataAnalysis.py"
%run "/Users/johny/Documents/nEXO/fits.py"
#SET number of waveforms to be read, cutoff freq for filter, and range of pulse
N = 20000; cut = 0.03; phsmin = 500; phsmax = 570; 
rawPHS_69V, filtPHS_69V = filterSpectraBin(N, cut, phsmin, phsmax); 
del dataFile;

#THEN, for filtered waveforms, to find the range of 4 largest peaks 
figure(); yr68,xr68,_= hist(filtPHS_68V, 200,histtype = 'step',label = '3 MHz lowpass filter, Vbias = -68 V');legend();
figure(); yr67,xr67,_= hist(filtPHS_67V, 200,histtype = 'step',label = '3 MHz lowpass filter, Vbias = -67 V');legend();
figure(); yr66,xr66,_= hist(filtPHS_66V, 200,histtype = 'step',label = '3 MHz lowpass filter, Vbias = -66 V');legend();
figure(); yr65,xr65,_= hist(filtPHS_65V, 200,histtype = 'step',label = '3 MHz lowpass filter, Vbias = -65 V');legend();
figure(); yr64,xr64,_= hist(filtPHS_64V, 200,histtype = 'step',label = '3 MHz lowpass filter, Vbias = -64 V');legend();


#THEN, to produce the plot/binning with 4 peaks, and FIND OUT expected values (by sight)
figure(); yf67,xf67,_= hist(filtPHS_67V, 200, range = (-0.081,.0),histtype = 'step',label = '3 MHz lowpass filter, Vbias = -67 V');legend();
figure(); yf70,xf70,_= hist(filtPHS_70V, 200, range = (-0.11,.01),histtype = 'step',label = '3 MHz lowpass filter, Vbias = -70 V');legend();
figure(); yf71,xf71,_= hist(filtPHS_71V, 200, range = (-0.13,.01),histtype = 'step',label = '3 MHz lowpass filter, Vbias = -71 V');legend();

#Lastly, for each peak (expected=(mu,sigma,A,...)):
figure(); yf315,xf315,_=hist(filtPHS_315V, 100, range=(-0.136,-0.021), histtype = 'step',label = '3 MHz filter, Vbias = -31.5 V, LED on'); legend();  
expected=(-0.03, 0.01, 1420, -0.06, 0.01, 520, -0.09, 0.01, 270, -0.12, 0.01, 140);
params315 = fourFitter(xf315,yf315,expected);
print('resolution and peak to peak for 31.5 V is: ');
resolution(params315);

#Trifitter example 
figure(); yn305,xn305,_=hist(filtPHS_noise_305V, 100, range=(-0.047,0), histtype = 'step',label = '3 MHz filter, Vbias = -30.5 V, LED off');legend();  
expected=(-0.008, 0.01, 820, -0.02, 0.01, 190, -0.04, 0.01, 16);
params305 = triFitter(xn305,yn305,expected);
print('resolution and peak to peak for 30.5 V is: ');
resolution_3peak(params305);

Routine for single peak fits: 
"""

# **For linear fits** 
#m,b = np.polyfit(Vbias,deltas,1)
#yList = []
#for i in range(len(Vbias)):
    #yFit = m*Vbias[i]+b
    #yList.append(yFit)
#plt.figure(); plt.scatter(Vbias, deltas); plt.plot(Vbias, yList, linestyle='solid')
#print('Vbreak is: '); 
#print -b/m + 1