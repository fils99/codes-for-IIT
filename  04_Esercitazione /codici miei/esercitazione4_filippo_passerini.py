#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 16:42:18 2022

@author: filippopasserini
"""

import time

start = time.time()


from funzioni4_filippo_passerini import dft, verifica_dft, idft, verifica_idft

import numpy as np
import matplotlib.pylab as plt
import math as math
from scipy import fftpack as f


print("\n Esercizi 1, 2, 3 \n")

v1=np.array([2,3])
v2=np.array([2,4,5])
v3=np.array([3,4,5,6])

V1=f.fft(v1)
V2=f.fft(v2)
V3=f.fft(v3)

W1=dft(v1)
W2=dft(v2)
W3=dft(v3)

verifica1=verifica_dft(V1,W1)
verifica2=verifica_dft(V2,W2)
verifica3=verifica_dft(V3,W3)

w1=idft(V1)
w2=idft(V2)
w3=idft(V3)

verifica_i1=verifica_idft(v1,w1)
verifica_i2=verifica_idft(v2,w2)
verifica_i3=verifica_idft(v3,w3)


print("\n Esercizio 4 \n")

X=np.arange(0,5)
x=idft(X)
y=np.zeros(len(x),dtype=complex)
for n in range(0,len(x)):
    y[n]=np.exp(4j*math.pi*n/5)*x[n]
Y=dft(y)
Y_prev=np.array([3,4,0,1,2])  # trasf prevista, per la proprietà di traslazione in frequenza
verifica_Y=verifica_dft(Y,Y_prev) 


print("\n Esercizio 5 \n")

t=np.linspace(0,1,50)
x1=np.zeros([50,3])
x2=np.zeros([50,3])
x3=np.zeros(50)
x4=np.zeros(50)
for i in range(0,50):
    for k in range(0,3):
        x1[i,k]=2018*np.sin(2*(k+1)*math.pi*t[i]) # x1 e x2 sono due matrici; colonna 0 corrisponde a k=1 ecc
        x2[i,k]=np.cos(2*(k+1)*math.pi*t[i])
        x3[i]=np.cos(2*math.pi*t[i])+np.cos(6*math.pi*t[i]+np.cos(10*math.pi*t[i]))
        x4[i]=np.exp(2*math.pi*1j*t[i])
                  
x11=x1[:,0]
x12=x1[:,1]
x13=x1[:,2]
x21=x2[:,0]
x22=x2[:,1]
x23=x2[:,2]

X11=dft(x11)
X12=dft(x12)
X13=dft(x13)
X21=dft(x21)
X22=dft(x22)
X23=dft(x23)
X3=dft(x3)
X4=dft(x4)

X11_ass=np.abs(X11)
X12_ass=np.abs(X12)
X13_ass=np.abs(X13)
X21_ass=np.abs(X21)
X22_ass=np.abs(X22)
X23_ass=np.abs(X23)
X3_ass=np.abs(X3)
X4_ass=np.abs(X4)

plt.stem(t,X11_ass)
plt.show()
plt.plot(t,X12_ass)
plt.show()
plt.plot(t,X13_ass)
plt.show()
plt.plot(t,X21_ass)
plt.show()
plt.plot(t,X22_ass)
plt.show()
plt.plot(t,X23_ass)
plt.show()
plt.plot(t,X3_ass)
plt.show()
plt.plot(t,X4_ass)
plt.show()

# si nota la linearità della trasformata di Fourier

print("\n Esercizio 6 \n")

s1=np.zeros(50)
s2=np.zeros(50)
s1[1]=10
s2[1]=10
s2[49]=10

S1=dft(s1)
S2=dft(s2)
S1_ass=np.abs(S1)
S2_ass=np.abs(S2)


plt.plot(t,S1_ass)
plt.show()
plt.plot(t,S2_ass)
plt.show()

# verifica identità parseval       
norm_quad_s1=np.linalg.norm(s1)**2
norm_quad_s2=np.linalg.norm(s2)**2
norm_quad_S1=np.linalg.norm(S1)**2
norm_quad_S2=np.linalg.norm(S2)**2
rap1=norm_quad_S1/norm_quad_s1
rap2=norm_quad_S2/norm_quad_s2


# si nota inoltre che la prima trasformata è cost, perchè il segnale ha una sola componente
# il secondo segnale invece ha due componenti, e dunque può essere visto come un segnale periodico con un'unica frequenza


print("Il tempo utilizzato per eseguire il codice è il seguente")

end = time.time()

print(end - start)

