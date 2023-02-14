#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 16:25:44 2022

@author: filippopasserini
"""

import time

start = time.time()




import numpy as np
import matplotlib.pylab as plt
import math as math
import scipy.io as sio
y=np.array(plt.imread('cameraman.jpg'), dtype='float64') # con float gli dico che è un'immagine reale
plt.imshow(y)
plt.show()


# usando base canonica, avremo base 64x64
# ogni coefficiente ottenuto dal prodotto scalare è un pixel
# solo L=3096 coefficienti vengono ricevuti
#generare maschera che butta 1000 coefficienti


dim=y.shape
print('le dimensioni dell immagine sono',dim)


control=0
for k in range(0,64):
    for i in range(0,64):
        if y[k,i]==math.floor(y[k,i]):
            control+=0
        else:
            control+=1
if control==0:
    print('i coefficienti dell immagine sono interi')
else:
    print('i coefficienti dell immagine non sono interi')
    
control2=0
for k in range(0,64):
    for i in range(0,64):
        if 0<=y[k,i]<=255:
            control2+=0
        else:
            control2+=1
if control2==0:
    print('i coefficienti dell immagine sono compresi fra 0 e 255')
else:
    print('i coefficienti dell immagine non sono compresi fra 0 e 255')

c=np.reshape(y,(4096,),order='F')
# c=np.zeros(4096)   modo alternativo
# i=0
# for k in range(0,64):
#     for j in range(0,64):
#             c[i]=y[j,k].copy()
#             i+=1

coef=np.zeros(4096)
x=np.random.permutation(4096)
for k in range(0,3096):  # oppure qualunque numero al posto di 3096
    indice=x[k]       
    coef[indice]=c[indice].copy()
    
       
C=np.reshape(coef,(64,64),order='F')
# C=np.zeros((64,64)) modo alternativo
# p=0   
# for k in range(0,64):
#     for j in range(0,64):
#         C[j,k]=coef[p]
#         p+=1

plt.imshow(C)
plt.show()



err_random=np.linalg.norm(y-C)/np.linalg.norm(y)   
print('l errore relativo dopo eliminazione casuale è',err_random)


B=sio.loadmat('DCT.mat')['B']


for k in range(0,4):
    b=np.reshape(B[k, :, :], [64, 64])
    plt.imshow(b)
    plt.show()

for k in range(64,68):
    b=np.reshape(B[k, :, :], [64, 64])
    plt.imshow(b)
    plt.show()

sopravvis=np.zeros(4096)   
for k in range(0,2048):    
        sopravvis[k]=c[k].copy()
    

S=np.reshape(sopravvis,(64,64),order='F')

# S=np.zeros((64,64))   modo alternativo
# t=0   
# for k in range(0,64):
#     for j in range(0,64):
#         S[j,k]=sopravvis[t]
#         t+=1

plt.imshow(S)
plt.show()

    

err_metà=np.linalg.norm(y-S)/np.linalg.norm(y)   
print('l errore relativo dopo aver eliminato la prima metà di coefficienti è',err_metà)

c_DCT=np.zeros(4096)
for u in range(0,4096):
        b=np.reshape(B[u, :, :], [64, 64])
        c_DCT[u]=np.sum(y*b)
       


sopravvis_DCT=np.zeros(4096)   
for k in range(0,2048):
        sopravvis_DCT[k]=c_DCT[k].copy()


I_rec=np.zeros((64,64))

for w in range(0,4096):
        b=np.reshape(B[w, :, :], [64, 64])
        I_rec+=sopravvis_DCT[w]*b
        

plt.imshow(I_rec)
plt.show()

err_metà_DCT=np.linalg.norm(y-I_rec)/np.linalg.norm(y)   
print('l errore relativo dopo aver eliminato la prima metà di coefficienti DCT è',err_metà_DCT)

if err_metà<err_metà_DCT:
    print('si ottiene un"approx migliore con la base canonica')
elif err_metà>err_metà_DCT:
    print('si ottiene un"approx migliore con la base DCT')
else:
    print('con le due basi si ottiene la stessa approx')
        
c_DCT_ass_sort=sorted(np.abs(c_DCT))

c_DCT_ass_ord=np.zeros(4096)
e=4095
for k in range(0,4096):
    c_DCT_ass_ord[k]=c_DCT_ass_sort[e]
    e-=1

sopravvis_DCT_bis=np.zeros(4096)   
for k in range(0,2048):
        sopravvis_DCT_bis[k]=c_DCT_ass_ord[k].copy()

# per ricordarmi dei segni di c_DCT
m=0
marker=np.zeros(4096)
a=np.zeros(4096)
for j in range(0,4096):
      if c_DCT[m]==np.abs(c_DCT[m]):
          marker[j]=1
          m+=1
      elif c_DCT[m]==-np.abs(c_DCT[m]):
              marker[j]=-1
              m+=1



c_ricomp=np.zeros(4096)
n=0
for r in range(0,4096):
      for n in range(0,2048):
        if sopravvis_DCT_bis[n]==np.abs(c_DCT[r]):
            c_ricomp[r]=np.abs(c_DCT[r])*marker[r]           


I_rec2=np.reshape(c_ricomp,(64,64))  

I_rec2=np.zeros((64,64))
for w in range(0,4096):
          b=np.reshape(B[w, :, :], [64, 64])
          I_rec2+=c_ricomp[w]*b
        

plt.imshow(I_rec2)
plt.show()

err_metà_DCT_bis=np.linalg.norm(y-I_rec2)/np.linalg.norm(y)
print('l errore relativo dopo aver eliminato la metà dei coefficienti DCT di modulo minore è',err_metà_DCT_bis)

# ripeto quest'ultimo procedimento, ma con la base canonica


coef_sort=np.sort(np.abs(c))[::-1]
ind = np.argsort(np.abs(c))[::-1]  # mi dà indici corrispondenti al riarrangiamento
p=int(1000)
plt.plot(np.arange(p),coef_sort[:p])
plt.show()


I_approx4 = np.zeros((y.shape[0],y.shape[1]))
# for ii in range(v):
I_approx4 = np.reshape(coef_sort,(64,64))
    
diff = y - I_approx4
plt.imshow(I_approx4)
plt.show()
relative_error4 = np.sqrt(np.sum(diff*diff))/np.sqrt(np.sum(y*y))
print('Errore migliori k base canonica:'+ str(relative_error4))



print("Il tempo utilizzato per eseguire il codice è il seguente")

end = time.time()

print(end - start)

    