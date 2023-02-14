#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 16:56:52 2022

@author: filippopasserini
"""


from  scipy  import fftpack as f
import numpy as np
from scipy.linalg import circulant



def conv_circ(segnale1,segnale2):
    C=circulant(segnale1)
    conv=np.dot(C,segnale2)
    return conv
 
    
def conv_fft(segnale1,segnale2):
    X1=f.fft(segnale1)
    X2=f.fft(segnale2)
    conv=np.real(f.ifft(X1*X2))
    return conv
