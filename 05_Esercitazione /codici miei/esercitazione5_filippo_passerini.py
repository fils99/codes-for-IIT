#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 16:28:26 2022

@author: filippopasserini
"""



import numpy as np
import matplotlib.pylab as plt
from scipy import fftpack as f
import scipy.io as sio
import scipy.io.wavfile
from numpy import linalg as la
import math
import itertools


print('\n esercizio 1 \n')

x=np.squeeze(sio.loadmat('frequencyRepresentation.mat')['x'])

plt.plot(x[0:500]) # non si riesce a distinguere quante frequenze ci stanno
plt.xlabel('Tempo')
plt.ylabel('x')         
plt.show()



scipy.io.wavfile.write('sound1.wav',8000,x)

X=f.fft(x)
plt.plot(np.abs(X)) # vedo che ho due frequenze
plt.xlabel('Frequenze')
plt.ylabel('|X|')
plt.show()


print('\n esercizio 2 \n')
y=np.zeros(128)
y[0:64]=1
Y=f.fft(y)

Y0=Y.copy()
Y1=Y.copy()
Y2=Y.copy()
Y3=Y.copy()
Y4=Y.copy()
Y5=Y.copy()
Y6=Y.copy()
Y7=Y.copy()
Y8=Y.copy()
Y9=Y.copy()

# for r in range(0,10):
    
Y0[1:]=0
Y1[2:]=0
Y2[3:]=0
Y3[4:]=0
Y4[5:]=0
Y5[6:]=0
Y6[7:]=0
Y7[8:]=0
Y8[9:]=0
Y9[10:]=0

y0=f.ifft(Y0)
y1=f.ifft(Y1)
y2=f.ifft(Y2)
y3=f.ifft(Y3)
y4=f.ifft(Y4)
y5=f.ifft(Y5)
y6=f.ifft(Y6)
y7=f.ifft(Y7)
y8=f.ifft(Y8)
y9=f.ifft(Y9)

plt.plot(y0)
plt.plot(y1)
plt.plot(y2)
plt.plot(y3)
plt.plot(y4)
plt.plot(y5)
plt.plot(y6)
plt.plot(y7)
plt.plot(y8)
plt.plot(y9)  # quello con più oscillazioni
plt.xlabel('Tempo')
plt.ylabel('|y_i|')
plt.show()

print('\n esercizio 3 \n')
N=1000
t=np.arange(0,N,1)
z1=np.sin(80*math.pi*t/N)
z2=np.sin(160*math.pi*t/N)
z3=np.sin(320*math.pi*t/N)     

zeri=np.zeros(1000)   

z_prov=list(itertools.chain(z1,zeri,z2,zeri,z3))
z=np.array(z_prov)   # lo trasformo da lista ad array
scipy.io.wavfile.write('sound2.wav',5000,z) # z in pratica è la somma di 3 segnali

Z=f.fft(z)
plt.plot(np.abs(Z))  # vediamo la dispersione spettrale
plt.xlabel('Frequenze')
plt.ylabel('|Z|')
plt.show()

Zapprox=Z.copy()
Zapprox[600:4399]=0
plt.plot(np.abs(Zapprox))
plt.xlabel('Frequenze')
plt.ylabel('|Zapprox|')
plt.show

zapprox=f.ifft(Zapprox)
zapprox2=np.array(zapprox,dtype=float) # lo converto in reale per poterlo salvare come file audio
scipy.io.wavfile.write('sound3.wav',5000,zapprox2)

# la differenza è dovuta al fatto che in pratica abbiamo perso una frequenza

err_real=np.max((zapprox-z).real)
err_imag=np.max((zapprox-z).imag)
err_relativo=la.norm(zapprox-z)/la.norm(z)

print('\n la massima differenza fra le parti reali è   \n',err_real)
print('\n la massima differenza fra le parti immaginarie è   \n',err_imag)
print('\n l"errore relativo è   \n',err_relativo)






