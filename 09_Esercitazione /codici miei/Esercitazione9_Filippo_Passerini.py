#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 16:20:28 2022

@author: filippopasserini
"""
import time
from  scipy  import fftpack as f
import numpy as np
import matplotlib.pylab as plt
from numpy import mean

start = time.time()
#Esercizio 1
print('\n Esercizio 1 \n')

x=np.array(plt.imread('barbara.png'),dtype='float64')
plt.imshow(x,cmap='gray')
plt.show()

def VisFourier(immagine,k):
    trasf=f.fft2(immagine)
    trasf_centrata=f.fftshift(trasf)
    log=np.log(np.abs(trasf_centrata)+1) # +1, perchè senno verrebbero fuori coeff negativi]]]
    plt.imshow(np.abs(trasf_centrata))
    if k==0:
        plt.title('no filtraggio')
    elif k==1:
        plt.title('post filtro passa basso')
    elif k==2:
        plt.title('post filtro passa alto')   
    elif k==3:
        plt.title('post filtraggio energia media') 
    plt.show()
    plt.imshow(log)
    if k==0:
        plt.title('no filtraggio')
    elif k==1:
        plt.title('post filtro passa basso')
    elif k==2:
        plt.title('post filtro passa alto')
    elif k==3:
        plt.title('post filtraggio energia media')
    plt.show()
    return trasf
# coef di maggiore intensità stanno al centro (frequenze basse), tipico delle immagini naturali, che sono piuttosto regolari
# con filtro passa basso quindi, l'immagine sarà circa preservata

k=0 # flag per i plot
X=VisFourier(x,k)


#Esercizio 2
print('\n Esercizio 2 \n')

k=1
def filtro_LP(immagine):
    N=int(immagine.shape[0])
    n=int(N/4)
    trasf=f.fftshift(f.fft2(immagine))
    trasf_mod=np.zeros((N,N))
    trasf_mod[n:int(3*n),n:int(3*n)]=np.ones((int(N/2),int(N/2)))
    trasf_mod=np.multiply(trasf_mod,trasf)
    im_filt=f.ifft2(f.ifftshift(trasf_mod))
    plt.imshow(np.real(im_filt),cmap='gray')
    plt.title('post filtro passa basso, 50% frequenze')
    plt.show()
    return im_filt

filtr=filtro_LP(x)
VisFourier(filtr,k)

def filtro_LP_alfa(immagine,alfa):
    N=int(immagine.shape[0])
    n=N*(1-alfa)/2
    n_arr=round(n,0)
    m=N-n
    m_arr=round(m,0)
    M=int(np.ceil(N*alfa))
    trasf=f.fftshift(f.fft2(immagine))
    trasf_mod=np.zeros((N,N))
    delta=M-(m_arr-n_arr)
    if delta!=0: # per evitare che dimensioni delle matrici non coincidano a causa di errori numerici 
        M=int(M-delta) 
    else:
        M=M
    trasf_mod[int(n_arr):int(m_arr),int(n_arr):int(m_arr)]=np.ones((M,M))
    trasf_mod=np.multiply(trasf_mod,trasf)
    im_filt=f.ifft2(f.ifftshift(trasf_mod))
    plt.imshow(np.real(im_filt),cmap='gray')
    plt.title('post filtro passa basso')
    plt.show()
    return im_filt

alfa_vec=np.array([0,0.05,0.1,0.2,0.7])

for ii in range(alfa_vec.shape[0]):
    alfa=alfa_vec[ii]
    filtr_2=filtro_LP_alfa(x,alfa)
    VisFourier(filtr_2,k)



#Esercizio 3
print('\n Esercizio 3 \n')


def filtro_HP_alfa(immagine,alfa):
    N=int(immagine.shape[0])
    n=N*(1-alfa)/2
    n_arr=round(n,0)
    m=N-n
    m_arr=round(m,0)
    M=int(np.ceil(N*alfa))
    trasf=f.fftshift(f.fft2(immagine))
    trasf_mod=np.ones((N,N))
    delta=M-(m_arr-n_arr)
    if delta!=0:
        M=int(M-delta)
    else:
        M=M
    trasf_mod[int(n_arr):int(m_arr),int(n_arr):int(m_arr)]=np.zeros((M,M))
    trasf_mod=np.multiply(trasf_mod,trasf)
    im_filt=f.ifft2(f.ifftshift(trasf_mod))
    plt.imshow(np.real(im_filt),cmap='gray')
    plt.title('post filtro passa alto')
    plt.show()
    return im_filt

k=2
for jj in range(alfa_vec.shape[0]):
    alfa=alfa_vec[jj]
    filtr_3=filtro_HP_alfa(x,alfa)
    VisFourier(filtr_3,k)

k=3
c=5 # inserire valore qualunque
media=mean(np.abs(X))
soglia=media*c
X_en=np.zeros((X.shape[0],X.shape[1]),dtype=complex)
for i in range(0,X.shape[0]):
    for j in range(0,X.shape[1]):
        if np.abs(X[i,j])>soglia:
            X_en[i,j]=X[i,j]

filtr_3=f.ifft2(X_en)           
VisFourier(filtr_3,k)       
plt.imshow(np.real(filtr_3),cmap='gray')
plt.title('post filtraggio energia media')
plt.show()
end=time.time()

print('il tempo (in s) impiegato dal codice è ',end-start)


