#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 12:14:27 2022

@author: filippopasserini
"""

import time
from  scipy  import fftpack as f
import numpy as np
import matplotlib.pylab as plt
from numpy import mean
import math as math
import scipy.io.wavfile
import scipy.io as sio

start = time.time()
#Esercizio 1
print('\n Esercizio 1 \n')


t=np.linspace(-25,25,1000)
x=(np.sinc(t))**2
plt.plot(t,x)
plt.xlabel('samples')
plt.ylabel('sinc')
plt.show()

X=f.fft(x)
normFrequ = np.arange(1,X.size+1,dtype=float)/float(X.size)
plt.plot(normFrequ,np.abs(X)) # è filtro passa basso 
plt.xlabel('normFrequ')
plt.ylabel('magnitude dft di sinc')
plt.show()

def campionamento(segnale,pc): # pc = passo di campionamento
    # segn_camp=np.zeros(int(len(segnale)/pc))
    segn_camp=segnale[0:len(segnale):pc]
    trasf_camp=f.fft(segn_camp)
    normFrequ = np.arange(1,trasf_camp.size+1,dtype=float)/float(trasf_camp.size)
    plt.plot(normFrequ,np.abs(trasf_camp))
    plt.xlabel('normFrequ')
    plt.ylabel('magnitude dft campionata')
    plt.show()
    return [segn_camp,trasf_camp]

def ricostruzione(trasf_camp,segnale_orig):
    L=int(len(trasf_camp))
    L2=int(L/2)
    L3=int(len(segnale_orig))
    M=L3-L   
    fattore_scal=float(len(f.fft(segnale_orig))/len(trasf_camp))
    z=np.zeros((M))
    trasf_zero=np.concatenate((trasf_camp[0:L2],z,trasf_camp[L2:L]))
    if M!=0:
       segn_rico=np.real(f.ifft(trasf_zero))*fattore_scal
    else:
       segn_rico=np.real(f.ifft(trasf_camp))*fattore_scal
    plt.plot(segn_rico)
    plt.xlabel('samples')
    plt.ylabel('segnale ricostruito')
    plt.show()
    return segn_rico
    

[x1,X1]=campionamento(x,5)
[x2,X2]=campionamento(x,10)
[x3,X3]=campionamento(x,20)

x1_rico=ricostruzione(X1,x)
x2_rico=ricostruzione(X2,x)
x3_rico=ricostruzione(X3,x)
         

# non riesco a recuperare segnale nel terzo caso; x3 non è a banda limitata

#Esercizio 2
print('\n Esercizio 2 \n')

t=np.linspace(0,999,1000)

y1=np.sin(80*math.pi*t/1000)
y2=np.sin(160*math.pi*t/1000)
y3=np.sin(320*math.pi*t/1000)
z=np.zeros(1000)
y=np.concatenate((y1,z,y2,z,y3))
plt.xlabel('samples')
plt.ylabel('segnale y')
plt.plot(y)
plt.show()

scipy.io.wavfile.write('segnaley.wav',9000,y)
Y=f.fft(y)
normFrequ = np.arange(1,Y.size+1,dtype=float)/float(Y.size)
plt.plot(normFrequ,np.abs(Y))
plt.xlabel('normFrequ')
plt.ylabel('magnitude dft di y')
plt.show()


[y1,Y1]=campionamento(y,2)
y1_rico=ricostruzione(Y1,y)
scipy.io.wavfile.write('segnaley_rico1.wav',9000,y1_rico)
# segnale riprodotto bene, ma il volume è un po' più basso
# frequenze leggermente aumentate

[y2,Y2]=campionamento(y,10)
y2_rico=ricostruzione(Y2,y)
scipy.io.wavfile.write('segnaley_rico2.wav',9000,y2_rico)
# segnale riprodotto male
# effetto di aliasing; non sono soddisfatte ipotesi del th shannon, si perde info


#Esercizio 3
print('\n Esercizio 3 \n')


j=np.squeeze(sio.loadmat('jingle.mat')['jingle'])
scipy.io.wavfile.write('jingle.wav',44000,j)
plt.plot(j)
plt.xlabel('samples')
plt.ylabel('jingle')
plt.show()

J=f.fft(j)
normFrequ = np.arange(1,J.size+1,dtype=float)/float(J.size)
plt.plot(normFrequ,np.abs(J))
plt.xlabel('normFrequ')
plt.ylabel('magnitude dft jingle')
plt.show()

[j1,J1]=campionamento(j,2)
j1_rico=ricostruzione(J1,j)
scipy.io.wavfile.write('jingle_rico1.wav',44000,j1_rico)
# si sente bene

[j2,J2]=campionamento(j,10)
j2_rico=ricostruzione(J2,j)
scipy.io.wavfile.write('jingle_rico2.wav',44000,j2_rico)
# si sente male, il volume è molto più basso
# non sono soddisfatte ipotesi th shannon

t=np.linspace(-50,50,101)
h=math.pi/10*np.sinc(t/10)

# H=f.fft(h) altro modo
# prod=H*J
# j_filtrato=np.real((f.ifft(prod)))
j_filtrato=np.convolve(j,h,'valid')

[j3,J3]=campionamento(j_filtrato,10)
j3_rico=ricostruzione(J3,j_filtrato)
scipy.io.wavfile.write('jingle_rico3.wav',44000,j3_rico)


end=time.time()

print('il tempo (in s) impiegato dal codice è ',end-start)
    
    
    
    

