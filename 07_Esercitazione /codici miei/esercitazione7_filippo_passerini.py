#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 16:38:03 2022

@author: filippopasserini
"""


from  scipy  import fftpack as f
import numpy as np
import matplotlib.pylab as plt
import math
from funzioni7_filippo_passerini import conv_circ, conv_fft
from time import process_time



#Esercizio 1
print('\n Esercizio 1 \n')

a=np.array([1,2,3,0,0,0])
b=np.array([2,3,4,0,0,0])
conv1=conv_circ(a,b)
conv2=conv_fft(a,b)


#Esercizio 2
print('\n Esercizio 2 \n')

a_bis=np.linspace(0,1023,1024)
b_bis=np.linspace(1,1024,1024)

t1_start = process_time() 
conv1_bis=conv_circ(a_bis,b_bis)
t1_stop = process_time()
tempo1=t1_stop-t1_start

t2_start = process_time() 
conv2_bis=conv_fft(a_bis,b_bis)
t2_stop = process_time()
tempo2=t2_stop-t2_start

if tempo1<tempo2:
    print('\n la convoluzione circolare calcolata attraverso la moltiplicazione matriciale è più rapida rispetto a quella calcolata con la DFT \n')
elif tempo1>tempo2:
              print('\n la convoluzione circolare calcolata attraverso la DFT è più rapida rispetto a quella calcolata con la moltiplicazione matriciale \n')
else:
                  print('\n i due metodi utilizzati per calcolare la convoluzione circolare impiegano lo stesso tempo \n')


#Esercizio 3
print('\n Esercizio 3 \n')
t1=np.linspace(-50,50,101)
h=40*math.pi/1000*np.sinc(40*math.pi*t1/1000)
H=f.fft(h)
normFrequ = np.arange(1,H.size+1,dtype=float)/float(H.size)
H1=np.fft.fftshift(H)
plt.plot(normFrequ,np.abs(H1))
plt.show()

t=np.linspace(0,999,1000)
x1=np.sin(40*math.pi*t/1000)
x2=np.sin(320*math.pi*t/1000)



X1=f.fft(x1)
X2=f.fft(x2)
normFrequ=np.arange(1,X1.size+1,dtype=float)/float(X1.size)
X1_centrato=np.fft.fftshift(X1)
X2_centrato=np.fft.fftshift(X2)
plt.plot(normFrequ,np.abs(X1_centrato))
plt.plot(normFrequ,np.abs(X2_centrato))
plt.show() 
      
y1=np.convolve(h,x1,'valid')
y2=np.convolve(h,x2,'valid')
Y1=f.fft(y1)
Y2=f.fft(y2)
normFrequ=np.arange(1,Y1.size+1,dtype=float)/float(Y1.size)
Y1_centrato=np.fft.fftshift(Y1)
Y2_centrato=np.fft.fftshift(Y2)
plt.plot(normFrequ,np.abs(Y1_centrato))
plt.show()       
plt.plot(normFrequ,np.abs(Y2_centrato))
plt.show()

x=2*x1-10*x2
y_tot=np.convolve(h,x,'valid')
y_tot_bis=2*y1-10*y2

if (np.max(y_tot-y_tot_bis))<10**(-10):
    print('\n la linearità è verificata \n')
    
else:
    print('\n la linearità NON è verificata \n') 
    
    
x1d=np.roll(x1,100) # con roll fai rotolare il vettore, invece la prof ha proprio eliminato le prime 100 componenti
y1d=np.roll(y1,100)
y1d_bis=np.convolve(h,x1d,'valid')

if (np.max(y1d-y1d_bis))<10**(-10):
    print('\n la tempo-invarianza è verificata \n')
    
else:
    print('\n la tempo-invarianza NON è verificata \n') 
    

# quando fai prodotto fra x2 e h, poichè x2 ha freq non inclusa nel supporto della di h, vengono soppresse
    