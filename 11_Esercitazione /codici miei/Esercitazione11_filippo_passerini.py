#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 16:55:21 2022

@author: filippopasserini
"""

import time
import numpy as np
import matplotlib.pylab as plt
import scipy.signal as sissi

start = time.time()
#Esercizio 1
print('\n Esercizio 1 \n')


def downsampling(immagine,coef): # coef è il fattore di downsampling
    D=np.zeros((int(immagine.shape[0]/coef),int(immagine.shape[0])))
    for j in range(0,D.shape[0]):
        D[j,coef*j]=1
    D_bis=D.transpose()
    down_prov=np.dot(D,immagine) # provvisorio
    down=np.dot(down_prov,D_bis)
    plt.imshow(down,cmap='gray')
    plt.title('immagine post downsampling')
    plt.show()
    return down

def sovracampionamento4(immagine):  
    t=0
    S=np.zeros((int(immagine.shape[0]*2),int(immagine.shape[0])))
    while t<S.shape[1]:    
        for i in range(0,int(S.shape[0]),2): 
                    S[i,int(t)]=1
                    S[i+1,int(t)]=1
                    t+=1      
    S_bis=S.transpose()
    sovra_prov=np.dot(S,immagine) # provvisorio
    sovra=np.dot(sovra_prov,S_bis)
    plt.imshow(sovra,cmap='gray')
    plt.title('immagine post downsampling e sovracampionamento')
    plt.show()
    return sovra

x=np.array(plt.imread('peppers.png'), dtype='float64') 
plt.imshow(x,cmap='gray')
plt.title('immagine originale')
plt.show()
F=1/25*np.ones((5,5))
x_fil=sissi.convolve2d(x,F,'same') # con same, all'uscita ho immagine di stesse dimensioni dell'originale
plt.imshow(x_fil,cmap='gray')
plt.title('immagine post filtro passa basso')        
plt.show()

N=np.random.normal(0,0.05,x.shape)
xn=x+N
plt.imshow(xn,cmap='gray')
plt.title('immagine sporcata')
plt.show()
x_denoise=sissi.convolve2d(xn,F,'same')
plt.imshow(x_denoise,cmap='gray')
plt.title('immagine ripulita')
plt.show()



#Esercizio 2
print('\n Esercizio 2 \n')


F2=np.array([[1/2],[0],[-1/2]])
x_fil_2=sissi.convolve2d(x,F2,'same')
plt.imshow(x_fil_2,cmap='gray') # mette in risalto i bordi orizzontali
plt.title('immagine post filtro passa alto')
plt.show()

F2bis=F2.transpose()
x_fil_2bis=sissi.convolve2d(x,F2bis,'same')
plt.imshow(x_fil_2bis,cmap='gray') # mette in risalto i bordi verticali
plt.title('immagine post filtro passa alto del secondo tipo')
plt.show()

grad=(x_fil_2**2+x_fil_2bis**2)**0.5
plt.imshow(grad,cmap='gray')
plt.title('grad') # dunque è un modo per estrarre i contorni
plt.show()

x_D2=downsampling(x,2)
x_fil_D2=downsampling(x_fil,2)

x_S4_D2=sovracampionamento4(x_D2)
x_fil_S4_D2=sovracampionamento4(x_fil_D2)


#Extra
print('\n Extra \n')

K1=np.array([[0,1,0],[1,-4,1],[0,1,0]])
x_fil3=sissi.convolve2d(x,K1,'same') 
plt.imshow(x_fil3,cmap='gray')
plt.title('immagine post filtro K1')        
plt.show()

K2=np.array([[2,1,0],[1,1,-1],[0,-1,-2]])
x_fil4=sissi.convolve2d(x,K2,'same') 
plt.imshow(x_fil4,cmap='gray')
plt.title('immagine post filtro K2')        
plt.show()

K3=np.array([[-2,-1,0],[-1,-1,1],[0,1,2]])
x_fil5=sissi.convolve2d(x,K3,'same') 
plt.imshow(x_fil5,cmap='gray')
plt.title('immagine post filtro K3')        
plt.show()



end=time.time()


print('il tempo (in s) impiegato dal codice è ',end-start)
    
