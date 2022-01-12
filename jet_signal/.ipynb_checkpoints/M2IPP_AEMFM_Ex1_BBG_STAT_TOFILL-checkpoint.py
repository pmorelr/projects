#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 20:11:34 2017

@author: P.Chary from R. Monchaux matlab scripts

Romain Monchaux
M2 IPP - Advanced Experimental Methods in FLuid Mechanics
Part 1: Average, standard deviation and probability density function 
        of a Gaussian white noise

"""

# Import required modules
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import scipy.signal
import scipy.integrate as spi


# Load signal
temp1 = scipy.io.loadmat('../signaux/sig1.mat')
t1 = temp1['t'].flatten()      # get time vector t1
Sig1 = temp1['Sig1'].flatten() # get signal vector
Fs1 = int(temp1['Fs'])         # get sampling frequency

## Question 1: plot raw signal:
plt.figure()
plt.plot(t1,Sig1,'.-')
plt.xlabel('temps (s)')
plt.ylabel('Signal')
plt.title('Sig1')
plt.show()


# Zoom in at different scales
xmin = 
xmax = 
plt.xlim([xmin,xmax/1000])
plt.savefig('../figures/BBG_de_tps_'+str(xmax)+'.png',  bbox_inches = 'tight')


## Question 2: Average and standar deviation estimations 
#              with a varying number of data points
NbPt = 2 ** np.array([XX,XX,XX])
average, standardev = [],[]
for i in range(len(NbPt)):
    Data=Sig1[:NbPt[i]]
    average.append(np.mean(Data))
    standardev.append(np.std(Data))

plt.figure()
plt.plot(NbPt,average,'-or')
plt.plot(NbPt,standardev,'-sb')
plt.plot([2,2**20],[0,0 ],'-.k')
plt.plot([2,2**20],[1,1],'-.k')
plt.xlabel('Nombre de points')
plt.ylabel('Valeur du moment')
plt.title('Sig1')
plt.legend({'average','Standard deviation'})
plt.xscale('log')

plt.savefig('../figures/BBG_Moy_de_N.png',  bbox_inches = 'tight')


## Question 3: estimation de la densité de probabilités:
def Gauss(x):
 return 1 / (np.sqrt(2 * np.pi)) * np.exp(-(x ** 2) / 2)

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

# With a single NbPt/NbBin combination
NbPt = XX
NbBin = XX
Bin,PDF = Calc_PDF(Sig1[:NbPt],NbBin)
BG = np.linspace(-5,5,2**20)

fig, ax = plt.subplots(1, 3, figsize=(12, 4), sharey=False)
fig.tight_layout(h_pad=2)
ax[0].plot(Bin,PDF)
ax[0].plot(BG,Gauss(BG),'-.k')
ax[0].set_xlabel('Signal Values')
ax[0].set_ylabel('PDF')

ax[1].plot(Bin,PDF)
ax[1].plot(BG,Gauss(BG),'-.k')
ax[1].set_xlabel('Signal Values')
ax[1].set_ylim([0.25,0.45])
ax[1].set_xlim([-1,1])

ax[2].plot(Bin,PDF)
ax[2].plot(BG,Gauss(BG),'-.k')
ax[2].set_xlabel('Signal Values')
# ax[2].set_ylabel('PDF')
ax[2].set_yscale('log')

fig.savefig('../figures/BBG_PDF_0.png', bbox_inches = 'tight')


NbPt = 2**XX
NbBin = [XX, XX, XX, XX]

fig, ax = plt.subplots(2, 5, figsize=(12, 4), sharey=False)
fig.tight_layout(h_pad=2)
for i in range(0,len(NbBin)):
    Bin, PDF = Calc_PDF(Sig1[:NbPt], NbBin[i])
    ax[0,i%5].plot(Bin,PDF)
    ax[0,i%5].plot(BG,Gauss(BG),'-.k')
    ax[0,i%5].set_title('NbBin='+str(NbBin[i]))
    ax[0,i%5].set_xlabel('Signal Values')
    ax[0,i%5].set_ylabel('PDF')
    ax[0,i%5].set_ylim([0.25,0.45])
    ax[0,i%5].set_xlim([-1,1])
    
    ax[1,i%5].plot(Bin,PDF)
    ax[1,i%5].plot(BG,Gauss(BG),'-.k')
    ax[1,i%5].set_yscale('log')
    ax[1,i%5].set_xlabel('Signal Values')
    ax[1,i%5].set_ylabel('PDF')

fig.savefig('../figures/BBG_PDF.png',  bbox_inches = 'tight')


## Synthesis:
plt.figure()
plt.xlabel('Number of points')
plt.ylabel('Number of Bins')
plt.title('PDF convergence')
NbPt  = 2 ** np.array([10, 12, 16, 18, 20])
NbBin =    [25, 40, 65, 120, 200] 
plt.plot(NbPt, NbBin, 'ob',label = 'data points')
plt.plot(np.linspace(0, 10**6, 10**5), 2 * np.linspace(0, 10**6, 10**5)**(1/3),'-.k', label='y = 2 x^(1/3)')
plt.legend()
plt.savefig('../figures/BBG_synthese.png',  bbox_inches = 'tight')


## Spectrum:
NbPoint = 2**20 # number of points in the signal
Nfft = 2**XX # number of points used for the FFT calculation
Noverlap = Nfft / 2 # number of overlaping points
Nwindow = Nfft # window length
Fs = Fs1 # sampling frequency

F, Pxx = scipy.signal.welch(Sig1[:NbPoint], Fs, window = 'hanning', nperseg = Nwindow, noverlap = Noverlap, nfft = Nfft)

plt.figure()
plt.plot(F,Pxx)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density')
plt.title('L signal=2^{'+str(int(np.log(Nwindow)/np.log(2)))+'} - Nfft=2^{'+str(int(np.log(Nfft)/np.log(2)))+'}')
plt.text(100,2.5*10**(-3),'$\Delta f=$'+str(round(Fs/Nfft,3)))
plt.savefig('../figures/BBG_DSP_Nfft'+str(Nfft)+'.png',  bbox_inches = 'tight')


## Autocorrélation:
test0 = Sig1[:50000] # too shorten the calculation
test = test0 - np.mean(test0) # centering
testnorm = np.sum(test**2)
acor = np.correlate(test,test,mode="same") / testnorm
acor = acor[len(acor)//2:]
plt.figure()
plt.xlabel('$\Delta$ t')
plt.ylabel('Normalised autocorrelation')
delta_t=np.linspace(0,len(acor)/Fs1,len(acor))
plt.plot(delta_t,acor,'o-')
plt.xlim([-0.001,.01])
plt.savefig('../figures/BBG_Autocorr.png',  bbox_inches = 'tight')
