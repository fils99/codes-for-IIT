#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 16:26:08 2022

@author: filippopasserini
"""

from  scipy  import fftpack as f
import numpy as np
import matplotlib.pylab as plt
import scipy.io as sio
from scipy.ndimage.interpolation import shift


#Esercizio 1
print('\n Esercizio 1 \n')

x=np.squeeze(sio.loadmat('jingle.mat')['jingle'])
plt.plot(x)
plt.xlabel('tempo')
plt.ylabel('segnale')
plt.show()
X=f.fft(x)
normFrequ = np.arange(1,X.size+1,dtype=float)/float(X.size)
plt.plot(normFrequ,np.abs(X))
plt.xlabel('normFrequ')
plt.ylabel('|X|')
plt.show()

sio.wavfile.write('jingle.wav',44000,x)

n=np.random.normal(0,0.01,len(x))

xn=x+n
plt.plot(xn)
plt.xlabel('tempo')
plt.ylabel('segnale sporcato')
plt.show()
XN=f.fft(xn)
normFrequ = np.arange(1,XN.size+1,dtype=float)/float(XN.size)
plt.plot(normFrequ,np.abs(XN))
plt.xlabel('normFrequ')
plt.ylabel('|XN|')
plt.show()
sio.wavfile.write('jingle_sporcato.wav',44000,xn)


# ripulisco segnale con medie mobili
M=10
h=1/M*np.ones(M)
y=np.convolve(xn,h,'valid')
sio.wavfile.write('jingle_pulito.wav',44000,y)
Y=f.fft(y) 
normFrequ = np.arange(1,Y.size+1,dtype=float)/float(Y.size)
plt.plot(normFrequ,np.abs(Y))
plt.xlabel('normFrequ')
plt.ylabel('|Y|')
plt.show()


# ripulisco segnale con leaky integrator
lam=0.9
t=np.linspace(0,99,100)
u=np.ones(len(t))
h2=(1-lam)*lam**t*u
y2=np.convolve(xn,h2,'valid')
sio.wavfile.write('jingle_pulito_2.wav',44000,y2)
Y2=f.fft(y2) 
normFrequ = np.arange(1,Y2.size+1,dtype=float)/float(Y2.size)
plt.plot(normFrequ,np.abs(Y2))
plt.xlabel('normFrequ')
plt.ylabel('|Y2|')
plt.show()
# vedo che filtrando il segnale, levo delle frequenze in eccesso e riduco un po' l'ampiezza del segnale


#Esercizio 2
print('\n Esercizio 2 \n')

dates=np.squeeze(sio.loadmat('data.mat')['dates_ts'])
p=np.squeeze(sio.loadmat('data.mat')['price_ts'])

r=(p-shift(p,-1,cval=0))/p
plt.plot(r[:len(r)-1])
plt.xlabel('tempo')
plt.ylabel('ritorno dell"investimento')
plt.show()

# con medie mobili
M=50
h3=1/M*np.ones(M)
r_medio=np.convolve(r,h3,'valid')
r_medio_bis=np.convolve(r**2,h3,'valid')
sigma=((r_medio_bis-2*r_medio**2)/M+r_medio**2)**0.5
plt.plot(r_medio[:len(r_medio)-1])
plt.xlabel('tempo')
plt.ylabel('ritorno medio')
plt.show()
plt.plot(sigma[:len(sigma)-1])
plt.xlabel('tempo')
plt.ylabel('volatilità')
plt.show()

# con leaky integrator
lam2=0.94
t2=np.linspace(0,99,100)
u2=np.ones(len(t2))
h4=(1-lam2)*lam2**t2*u2
r_medio2=np.convolve(r,h4,'valid')
r_medio2_bis=np.convolve(r**2,h4,'valid')
sigma2=((r_medio2_bis-2*r_medio2**2)/M+r_medio2**2)**0.5
plt.plot(r_medio2[:len(r_medio2)-1])
plt.xlabel('tempo')
plt.ylabel('ritorno medio2')
plt.show()
plt.plot(sigma2[:len(sigma2)-1])
plt.xlabel('tempo')
plt.ylabel('volatilità2')
plt.show()

