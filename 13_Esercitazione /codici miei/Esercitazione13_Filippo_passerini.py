#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 16:32:16 2022

@author: filippopasserini
"""

import time
import numpy as np
import matplotlib.pylab as plt
import scipy.signal as sissi
from  scipy  import fftpack as f
import scipy.io as sio
import pywt as pywt
import scipy.io.wavfile


# decomposiz lineare e non lineare usando wavelet di 1 livello e 2 livelli

#Esercizio 1
print('\n Esercizio 1 \n')



x=np.array(plt.imread('barbara.png'),dtype='float64')
plt.imshow(x,cmap='gray')
plt.show()
cA, cD=pywt.dwt2(x,'haar')

plt.imshow(cA)
plt.title('approssimazione')
plt.show()

plt.imshow(cD[0])
plt.title('dettagli orizzontali') # vedo bordi orizzontali
plt.show()

plt.imshow(cD[1])
plt.title('dettagli verticali') # vedo bordi verticali
plt.show()

plt.imshow(cD[2])
plt.title('dettagli obliqui') # # vedo bordi diagonali
plt.show()

xrec=pywt.idwt2((cA, cD),'haar')

err_quad_medio=np.mean((x-xrec)**2)
print('errore quadratico medio:  ',err_quad_medio)

cD0_appr=np.zeros(cD[0].shape)
cD1_appr=np.zeros(cD[1].shape)
cD2_appr=np.zeros(cD[2].shape)

xrec1=pywt.idwt2((cA, (cD0_appr, cD1_appr, cD2_appr)),'haar') # butto gli obliqui
err_quad_medio1=np.mean((x*-xrec1)**2)
print('errore quadratico medio1:  ',err_quad_medio1)

xrec2=pywt.idwt2((cA, (cD[0], cD[1], cD2_appr)),'haar') # butto gli obliqui
err_quad_medio2=np.mean((x*-xrec2)**2)
print('errore quadratico medio2:  ',err_quad_medio2)

xrec3=pywt.idwt2((cA, (cD[0], cD1_appr, cD[2])),'haar') # butto i verticali
err_quad_medio3=np.mean((x-xrec3)**2)
print('errore quadratico medio3:  ',err_quad_medio3)

xrec4=pywt.idwt2((cA, (cD0_appr, cD[1], cD[2])),'haar') # butto gli orizzontali
err_quad_medio4=np.mean((x-xrec4)**2)
print('errore quadratico medio4:  ',err_quad_medio4)

# errore maggiore quando butto i dettagli verticali
# dettagli obliqui sono a bassa intensità, errore piccolo se li butto


# ora butto solo i coef meno espressi

soglia=np.mean(np.abs(cD))

indici_0=np.where(cD[0]>soglia)
indici_1=np.where(cD[1]>soglia)
indici_2=np.where(cD[2]>soglia)

cD0_appr[indici_0[0],indici_0[1]]=cD[0] [indici_0[0],indici_0[1]]
cD1_appr[indici_1[0],indici_1[1]]=cD[1] [indici_1[0],indici_1[1]]
cD2_appr[indici_2[0],indici_2[1]]=cD[2] [indici_2[0],indici_2[1]]

xrec5=pywt.idwt2((cA, (cD0_appr, cD1_appr, cD2_appr)),'haar') 
err_quad_medio5=np.mean((x-xrec5)**2)
print('errore quadratico medio5:  ',err_quad_medio5)

#Esercizio 6
print('\n Esercizio 6 \n')

cA2, cD1, cD2=pywt.wavedec2(x,'db4',level=2)

#Esercizio 7
print('\n Esercizio 7 \n')
xrec_db4=pywt.waverec2((cA2,cD2,cD1),'db4')
err_quad_medio_db4=np.mean((x-xrec_db4)**2) 

#Esercizio 8
print('\n Esercizio 8 \n')

media=np.mean(np.abs(cD2))

indici_d1=np.where(np.abs(cD1)>media)
indici_d2=np.where(np.abs(cD2)>media)
indici_a2=np.where(np.abs(cA2)>media)

cD1_appr=np.zeros(len(cD1))
cD1_appr[indici_d1[0],indici_d1[1]]=cD1[indici_d1[0],indici_d1[1]]

cD2_appr=np.zeros(len(cD2))
cD2_appr[indici_d2[0],indici_d2[1]]=cD2[indici_d2[0],indici_d2[1]]

cA2_appr=np.zeros(len(cA2))
cA2_appr[indici_a2[0],indici_a2[1]]=cD2[indici_a2[0],indici_a2[1]]

x_db4_approx=pywt.waverec2((cA2_appr,cD2_appr,cD1_appr),'db4')


           
#wavelet sono una buona base per le immagini naturali      