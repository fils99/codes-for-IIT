#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:42:17 2022

@author: filippopasserini
"""


import numpy as np
import matplotlib.pylab as plt
from scipy import fftpack as f
import scipy.io.wavfile
from funzioni4_filippo_passerini import dft
import math


#Esercizio 1
print('\nEsercizio 1\n')
x=np.zeros(128)
x[0:64] = 1

# plotto la mia dft
X_mia=dft(x)
normFrequ = np.arange(1,X_mia.size+1,dtype=float)/float(X_mia.size)
ass_mia=np.abs(X_mia)
plt.stem(normFrequ,ass_mia)
plt.xlabel('Normalized Frequencies')
plt.ylabel('modulo X_mia')
plt.show()
fase_mia=np.angle(X_mia)
plt.stem(normFrequ,fase_mia)
plt.xlabel('Normalized Frequencies')
plt.ylabel('fase X_mia')
plt.show()

# plotto la dft di python

X=f.fft(x)
ass=np.abs(X)
plt.stem(normFrequ,ass)
plt.xlabel('Normalized Frequencies')
plt.ylabel('modulo X')
plt.show()
fase=np.angle(X)
plt.stem(normFrequ,fase)
plt.xlabel('Normalized Frequencies')
plt.ylabel('fase X')
plt.show()



# quando la fase è vicina a un multiplo di pi, si introducono degli errori numerici
# conviene usare la DFT di python, perchè è più veloce ed evita errori numerici
# con il tuo codice, invece che zeri, ho quantità molto piccole ma non zero, e dunque si crea errore numerico
# errore numerico si vede nella fase, non nel modulo

#Esercizio 2
print('\nEsercizio 2\n')
t=np.linspace(0,3999,4000)
y=np.sin(80*math.pi*t/1000)+np.sin(100*math.pi*t/1000)
normFrequ = np.arange(1,y.size+1,dtype=float)/float(y.size)
Y=f.fft(y)
Y_ass=np.abs(Y)
plt.plot(normFrequ,Y_ass)
plt.xlabel('Normalized Frequencies')
plt.ylabel('modulo Y')
plt.show()

M=500
Y500=f.fft(y,M)
normFrequ = np.arange(1,M+1,dtype=float)/float(M)
Y_ass500=np.abs(Y500)
plt.plot(normFrequ,Y_ass500)
plt.xlabel('Normalized Frequencies')
plt.ylabel('modulo Y500')
plt.show()

M=50
Y50=f.fft(y,M)
normFrequ = np.arange(1,M+1,dtype=float)/float(M)
Y_ass50=np.abs(Y50)
plt.plot(normFrequ,Y_ass50)
plt.xlabel('Normalized Frequencies')
plt.ylabel('modulo Y50')
plt.show()


# se prendo multipli di 100 non perdo info rilevanti dello spettro, perchè le due freq sono 20 e 25, e 100 è il minimo comune multiplo
# allora usando meno campioni, ma comunque un numero multiplo 100, mantengo le info e ottimizzo, perchè diminuisco i calcoli

M=4000
y50=f.ifft(Y50)
Y50bis=f.fft(y50,M)
normFrequ = np.arange(1,M+1,dtype=float)/float(M)
Y_ass50bis=np.abs(Y50bis)
plt.plot(normFrequ,Y_ass50bis)
plt.xlabel('Normalized Frequencies')
plt.ylabel('modulo Y50bis')
plt.show()


# sopravvivono solo due picchi
# ricalcolando con 4000, vedo che ho perso info in precedenza; non riottengo lo spettro del segnale originario
# dunque con lo zero pudding non recupero informazione

#Esercizio 3
print('\nEsercizio 3\n')

C4 = 261.63
D4 = 293.66
F4 = 349.23
G4 = 392.00


t1=np.linspace(1,2399,2400)
t2=np.linspace(1,1599,1600)
t3=np.linspace(1,3199,3200)
t4=np.linspace(1,3999,4000)

x1=np.cos(2*math.pi*C4*t1/8000) # Do
x2=np.cos(2*math.pi*C4*t2/8000) # Do
x3=np.cos(2*math.pi*D4*t3/8000) # Re
x4=np.cos(2*math.pi*C4*t3/8000) # Do
x5=np.cos(2*math.pi*G4*t3/8000) # Sol
x6=np.cos(2*math.pi*F4*t4/8000) # Fa

M=500
normFrequ = np.arange(1,M+1,dtype=float)/float(M)
X1=f.fft(x1,n=M)
X1_ass=np.abs(X1)
plt.plot(normFrequ,X1_ass[0:500])
plt.xlabel('Frequencies')
plt.ylabel('modulo X1')
plt.show()

M=500
X2=f.fft(x2,n=M)
X2_ass=np.abs(X2)
plt.plot(normFrequ,X2_ass[:500])
plt.xlabel('Frequencies')
plt.ylabel('modulo X2')
plt.show()

M=500
X3=f.fft(x3,n=M)
X3_ass=np.abs(X3)
plt.plot(normFrequ,X3_ass[:500])
plt.xlabel('Frequencies')
plt.ylabel('modulo X3')
plt.show()

M=500
X4=f.fft(x4,n=M)
X4_ass=np.abs(X4)
plt.plot(normFrequ,X4_ass[:500])
plt.xlabel('Frequencies')
plt.ylabel('modulo X4')
plt.show()

M=500
X5=f.fft(x5,n=M)
X5_ass=np.abs(X5)
plt.plot(normFrequ,X5_ass[:500])
plt.xlabel('Frequencies')
plt.ylabel('modulo X5')
plt.show()

M=500
X6=f.fft(x6,n=M)
X6_ass=np.abs(X6)
plt.plot(normFrequ,X6_ass[:500])
plt.xlabel('Frequencies')
plt.ylabel('modulo X6')
plt.show()


z=np.concatenate((x1,x2,x3,x4,x5,x6))
Z=f.fft(z)
Z_ass=np.abs(Z)
lung=x1.size+x2.size+x3.size+x4.size+x5.size+x6.size
normFrequ = np.arange(1,lung+1,dtype=float)/float(lung)
plt.plot(normFrequ,Z_ass)
plt.xlabel('Frequencies')
plt.ylabel('modulo Z')
plt.show()

z_real=np.real(z)
scipy.io.wavfile.write('happy_birthday.wav',8000,z_real)
            

