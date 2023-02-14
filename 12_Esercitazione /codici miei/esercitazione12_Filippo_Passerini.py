#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 09:58:15 2022

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


start = time.time()

#Esercizio 1
print('\n Esercizio 1 \n')


x=np.zeros((512))
x[200]=50
X=f.fft(x)
plt.plot(f.fftshift(np.abs(X)))
plt.show()



y=np.zeros((512))
y[100]=50
Y=f.fft(y)
plt.plot(f.fftshift(np.abs(Y)))
plt.show()

# ottengo stessi risultati


X2=sissi.cwt(x,sissi.ricker,np.arange(1,48)) # trasformata a 48 livelli
# si può anche decidere il numero di livelli da analizzare
plt.imshow((X2))
plt.show()
Y2=sissi.cwt(y,sissi.ricker,np.arange(1,48)) 
plt.imshow((Y2))
plt.show()


z=np.zeros((512))
m=int(512/2) # scegliere valore fra 0 e len(z)/2
z[m]=100
z[int(len(z)-m)]=100
Z=sissi.cwt(z,sissi.ricker,np.arange(1,48)) 
plt.imshow((Z))
plt.show()



#Esercizio 2
print('\n Esercizio 2 \n')


seg=np.squeeze(sio.loadmat('chirp.mat')['x']) # squeeze per salvarlo come array 1D, estrae l'array dalla matrice
plt.plot(seg)
plt.show()


cA , cD=pywt.dwt(seg,'db4')
# in pratica ho fatto filtro passa basso e downsampling per avere cA
# filtro passa alto e downsampling per avere cD
A=pywt.upcoef('a',cA,'db4',take=seg.shape[0]) # take serve a imporre la lunghezza di A
D=pywt.upcoef('d',cD,'db4',take=seg.shape[0]) 
plt.plot(D)
plt.plot(A)
plt.show()
# upcoef serve a ricostruire il segnale singolarmente, o approssimazione o dettaglio
# upocef interpola info dei campioni vicini

seg_ric=pywt.idwt(cA,cD,'db4')
# idwt fa filtro inverso e upsampling; operazione ad albero al contrario
seg_ric2=A+D
err=np.max(seg_ric-seg_ric2)
print('l"errore max fra le due ricostruzioni è:  ',err)

cA3, cD3, cD2, cD1=pywt.wavedec(seg,'db4',level=3)
media=np.mean(np.abs(cD3))

indici_d1=np.where(np.abs(cD1)>media)
indici_d2=np.where(np.abs(cD2)>media)
indici_d3=np.where(np.abs(cD3)>media)
indici_a3=np.where(np.abs(cA3)>media)

cD1_appr=np.zeros(len(cD1))
cD1_appr[indici_d1]=cD1[indici_d1]

cD2_appr=np.zeros(len(cD2))
cD2_appr[indici_d2]=cD2[indici_d2]

cD3_appr=np.zeros(len(cD3))
cD3_appr[indici_d3]=cD3[indici_d3]

cA3_appr=np.zeros(len(cA3))
cA3_appr[indici_a3]=cA3[indici_a3]

seg_approx=pywt.waverec((cA3_appr,cD3_appr,cD2_appr,cD1_appr),'db4')
                   

err_quad_medio=np.mean((seg-seg_approx)**2)
print('l"errore quadratico medio dopo la compressione è  ', err_quad_medio)

scipy.io.wavfile.write('segnale.wav',9000,seg)
scipy.io.wavfile.write('segnale_approx.wav',9000,seg_approx)

end=time.time()


print('il tempo (in s) impiegato dal codice è ',end-start)
                       
