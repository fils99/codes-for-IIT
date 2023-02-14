#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 16:48:47 2022

@author: filippopasserini
"""


import numpy as np
import matplotlib.pylab as plt
import math as math
import scipy.io as sio
from scipy import fftpack as f

def dft(segnale):
    N=len(segnale)
    trasf=np.zeros(N,dtype=complex)
    for k in range(0,N):
        for n in range(0,N):
            trasf[k]+=complex(segnale[n])*np.exp(-2j*math.pi*k*n/N)
    return trasf

def verifica_dft(trasf1,trasf2):
    diff=trasf1-trasf2
    if np.max(abs(diff))<1e-10:
        print('la trasformata è la stessa')
    else:
        print('la trasformata NON è la stessa')
    return
    
    
def idft(trasf):
    N=len(trasf)
    segnale=np.zeros(N,dtype=complex)
    for n in range(0,N):
        for k in range(0,N):
            segnale[n]+=trasf[k]*np.exp(2j*math.pi*n*k/N)/N
    return segnale

def verifica_idft(segnale1,segnale2):
    diff=segnale1-segnale2
    if np.max(abs(diff))<1e-10:
        print('la anti-trasformata è la stessa')
    else:
        print('la anti-trasformata NON è la stessa')
    return


    
        
        
    
            
                  
          
    