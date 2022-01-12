#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 20:11:34 2017

@author: P.Chary from R. Monchaux matlab scripts

Romain Monchaux
M2 IPP - Advanced Experimental Methods in FLuid Mechanics
Part 2: Analysis of three su=ynthetic signals

"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import scipy.signal
import scipy.integrate as spi

# Load signals
temp2=scipy.io.loadmat('../signaux/sig2.mat')
temp3=scipy.io.loadmat('../signaux/sig3.mat')
temp4=scipy.io.loadmat('../signaux/sig4.mat')

################ Signal selection #################
number = 2 # set to 2, 3 or 4
#################################################

# Get data
temp_list=[_,_,temp2,temp3,temp4]
temp = temp_list[number]
# get time signal
t = temp['t'].flatten()
# get signal
Sig = temp['Sig' + str(number)].flatten()
# get sampling frequency
Fs = int(temp['Fs'])


## Plot time signals:
plt.figure()
plt.plot(t,Sig,'.-')
plt.xlabel('time (s)')
plt.ylabel('Signal')
plt.title(f'Sig{number}')
plt.show()

# Zoom on different signal parts:
xmin = XX
xmax = XX
plt.xlim([xmin,xmax/1000])
plt.savefig(f'../figures/Sig{number}_de_tps_{xmax}.png')

# Statistics:
average = np.mean(Sig)
standardev = np.std(Sig)

print(f'Signal {number} time average is {average} its standard deviation is {standardev}')

## Plot spectra:
# Number of points used in the signal
NbPoint = 2**XX
# Number of points used for FFT calculation
Nfft = 2**XX
# Number of points in the window
Nwindow = Nfft
# Number of overlapping points
Noverlap = Nfft/2

F, Pxx = scipy.signal.welch(Sig[:NbPoint],Fs,window='hanning',nperseg=Nwindow,noverlap=Noverlap,nfft=Nfft)

plt.figure()
plt.plot(F,Pxx)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density')
plt.title(f'DSP Sig{number} - Nfft = 2^{np.log(Nfft)/np.log(2)}')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1, Fs/2)

plt.savefig(f'../figures/Sig{number}_DSP_Nfft{Nfft}.png')

## Plot Probability Density Functions:
NbPt = 2**XX
NbBin = XX

def Gauss(x):
 return 1/(np.sqrt(2*np.pi))*np.exp(-(x**2)/2);


def Calc_PDF(SIG,NbBin):
    """ estimate the probability density function of the signal SIG 
    by integrating the histogram estimated over NbBin intervals
    The same can be done directly by using density option in np.histogram
    """
    # histogram estimation
    PDF, Bin = np.histogram(SIG, NbBin)
    Bin = Bin[:-1] + np.diff(Bin) / 2 # get bin centers
    # normalisize to probability density function
    Int = spi.trapz(PDF, Bin)
    PDF = PDF / Int
    return Bin, PDF

plt.figure()
BG = np.linspace(-5,5,2**20)
Bin, PDF = Calc_PDF(Sig[:NbPt],NbBin)
plt.subplot(2,1,1)
plt.plot(Bin,PDF)
plt.plot(BG,Gauss(BG),'-.k')
plt.xlabel('Signal Values')
plt.ylabel('PDF')
plt.title(f'PDF Sig{number}')
plt.ylim([10**(-6),10**0])
plt.xlim([-6,6])
plt.yscale('log')
plt.subplot(2,1,2)
plt.plot(Bin,PDF)
plt.plot(BG,Gauss(BG),'-.k')
plt.xlabel('Signal Values')
plt.ylabel('PDF')
plt.ylim([0,0.5])
plt.xlim([-2,2])

plt.savefig(f'../figures/Sig{number}_PDF.png')

## Plot autocorrellations:
test0 = Sig[:50000] #limitation sinon c'est trop long
test = test0 - np.mean(test0)
testnorm = np.sum(test**2)
acor = np.correlate(test,test,mode="same")/testnorm
acor = acor[len(acor)//2:]
plt.figure()
plt.xlabel('$\Delta$ t')
plt.ylabel('Normalised autocorrelation')
plt.title(f'Autocorr Sig{number}')
delta_t=np.linspace(0,len(acor)/Fs,len(acor))
plt.plot(delta_t,acor,'o-')
plt.xlim([-0.001,.01])
plt.savefig(f'../figures/Sig{number}_Autocorr.png')