#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 20:11:34 2017

@author: P.Chary from R. Monchaux matlab scripts

Romain Monchaux
M2 IPP - Advanced Experimental Methods in FLuid Mechanics
Part 3: real signal anaylisis

"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import scipy.io.wavfile
import scipy.signal


###########################################################################
# Question 1: 
# load and plot signal to test its stationarity
###########################################################################

# load signal
Fs,signal = scipy.io.wavfile.read('../signaux/gong_force.wav')
nbits = signal.dtype

Ts = 1/Fs # time step
Nsamples = len(signal) # sample number
time = np.linspace(0,Nsamples*Ts,Nsamples) # time vector
duration = round(Nsamples*Ts,1) # signal duration

# Split forcing/response
Fsig = signal[:,0]
Psig = signal[:,1]

# Plot time signals
plt.figure()
plt.subplot(2,1,1) # forcing
plt.plot(time,Fsig,'k-')
plt.title(f'Signal of duration {duration}s sampled at Fs={Fs/1000}kHz')

plt.xlabel('t (s)')
plt.ylabel('Forcing Amplitude (AU)')
#
plt.subplot(2,1,2) # pressure response
plt.plot(time,Psig,'k-')
plt.xlabel('t (s)')
plt.ylabel('pressure (Pa)')

plt.savefig('../figures/gong.png')


##########################################################################
# First statistical moments estimated over 1s 
###########################################################################
# sliding interval length (s)
Tint = XX
# Number of intervals
Nintervals = round(duration/Tint)
# time vector 
t1s = np.linspace(0,Nintervals*Tint,Nintervals)

# initialise averages and standard deviations
mean_f1s = np.zeros(Nintervals)
std_f1s = np.zeros(Nintervals)
mean_p1s = np.zeros(Nintervals)
std_p1s = np.zeros(Nintervals)
# number of samples per interval
Nsamples_int = int(Tint*Fs) 

# loop on intervals
for i in range(1,Nintervals): 
    mean_f1s[i] = np.mean(Fsig[(i-1)*Nsamples_int+1:i*Nsamples_int])
    std_f1s[i] = np.std(Fsig[(i-1)*Nsamples_int+1:i*Nsamples_int])
    mean_p1s[i] = np.mean(Psig[(i-1)*Nsamples_int+1:i*Nsamples_int])
    std_p1s[i] = np.std(Psig[(i-1)*Nsamples_int+1:i*Nsamples_int])

# plot
plt.figure()
plt.subplot(2,1,1)
plt.plot(t1s,mean_f1s,'ko-',label="mean forcing")
plt.plot(t1s,std_f1s,'r^-',label="standard deviation forcing")
plt.xlabel('t (s)')
plt.ylabel('forcing (UA)')
plt.title(f"Sliding moments over {Tint}s")
plt.legend(loc='best')
plt.grid()
#
plt.subplot(2,1,2)
plt.plot(t1s,mean_p1s,'ko-',label="mean pressure")
plt.plot(t1s,std_p1s,'r^-',label="standard deviation pressure")
plt.xlabel('t (s)')
plt.ylabel('pressure (Pa)')
plt.legend(loc='best')
plt.grid()
plt.savefig('../figures/gong_moy'+str(1000*Tint)+'.png')


###########################################################################
# Question 2: 
# power spectral density using pwelch
###########################################################################
  
# Full signal, no average
Nfft = 2**XX
# frequency resolution
df = Fs/Nfft
Nwindow = Nfft
Noverlap = int(Nfft/2)
# Spectrum estimation
F,Pxx = scipy.signal.welch(Psig[:Nwindow],Fs,window='hanning',nperseg=Nwindow,nfft=Nfft)

# plot
plt.figure()
plt.plot(F,Pxx)
plt.xlabel('Frequency (Hz)')
plt.ylabel('DSP (m$^2\cdot$s$^{-2}\cdot$Hz$^{-1}$)')
plt.grid()
plt.title('Power Spectral Density: gong')
plt.savefig(f'../figures/gong_freq_{Nfft}.png')

#plt.xlim([0,1000]) # zoom in frequency 10kHz (adjustable)
plt.yscale('log') # in log scale
plt.savefig('../figures/gong_spec_log.png')
plt.yscale('linear') # in linear scale
plt.savefig('../figures/gong_spec_lin.png')


##########################################################################
# Question 3: 
# Spectrogram using short time FFT
###########################################################################
# spectrogram of pressure response
# Frequency resolution
df = XX
# FFT point number
Nfft = int(Fs/df)
# window point number
Nwindow = Nfft
Noverlap = Nfft/2
# number of FFT segments
K = int(Nsamples/Nfft)
# time resolution
dt_specgram = Nwindow * Ts

# estimate spectrograms in each segment and corresponding times
SpecgramP = np.zeros((int(Nfft/2)+1,K)) # pressure
SpecgramF = np.zeros((int(Nfft/2)+1,K)) # Forcing

tspecgram = np.zeros(K)
# loop on segments
for i in range(1,K):
   # Pressure
   freq_specgramP,SpecgramP[:,i] = scipy.signal.welch(Psig[(i-1)*Nwindow+1:i*Nwindow],Fs,window='hanning',nperseg=Nwindow,noverlap=Noverlap,nfft=Nfft)
   # FOrcing
   freq_specgramF,SpecgramF[:,i] = scipy.signal.welch(Fsig[(i-1)*Nwindow+1:i*Nwindow],Fs,window='hanning',nperseg=Nwindow,noverlap=Noverlap,nfft=Nfft)
   # mid time of a segment
   tspecgram[i] = (i-0.5)*Nwindow*Ts 

# Plot for pressure
plt.figure()
plt.pcolor(tspecgram[1:],freq_specgramP[:-1]/1000,np.log(SpecgramP[:-1,1:])) #le [:,1:] permet d'enlever un terme nul qui donne un log valant -inf
plt.xlabel('t (s)')
plt.ylabel('f (kHz)')
plt.title('Pressure response spectrogram')
plt.colorbar()
# frequency zoom
plt.ylim([0,4])
plt.savefig('../figures/gong_spectrogram_P.png')

# plot for forcing 
plt.figure()
plt.pcolor(tspecgram[1:],freq_specgramF[:-1]/1000,np.log(SpecgramF[:-1,1:]))#le [:,1:] permet d'enlever un terme nul qui donne un log valant -inf
plt.xlabel('t (s)')
plt.ylabel('f (kHz)')
plt.title('Spectrogramme du signal de réponse en forçage')
plt.colorbar()
plt.ylim([0,4]) # zoom en fréquence
plt.savefig('../figures/gong_spectrogramme_F.png')
