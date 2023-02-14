#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 17:12:43 2022

@author: filippopasserini
"""


"""
Ex 1
"""
print('\n **** Esercizio  1 ****\n\n')

import numpy as np
import math
import matplotlib.pylab as plt
from funzioni2 import DFTinv, prodscal
import random as random

N=int(input('\n immettere valore di N   '))
C=np.zeros((N,N))
for n in range(0,N):
    for k in range(0,N):
        if k==0:
            C[n,k]=N**(-1/2)*math.cos((n+1/2)*k*math.pi/N)
            
        else:
            C[n,k]=(2/N)**(1/2)*math.cos((n+1/2)*k*math.pi/N)
print('\n matrice dei coefficienti C\n',C)          

a=np.round(np.dot(np.transpose(C),C),1)
print('\n prodotto fra C trasposto e C \n',a)

for k in range(0,N): # k ci dà la frequenza in pratica
    fig, ax=plt.subplots()
    plt.stem(C[:,k])

x=np.array([0, 0, 2, 3, 4, 0, 0, 0]).reshape((-1,1))
# oppure x=np.array([[0], [0], [2], [3], [4], [0], [0], [0]])

X=np.dot(np.transpose(C),x)
print('\n segnale trasformato \n',X)


x2=np.round(DFTinv(C,X,N),2)
print('\n segnale ricostruito \n',x2)


"""
Ex 2
"""
print('\n\n **** Esercizio  2 ****\n\n')


I1=np.floor(np.random.uniform(0,255,[64, 64]))
I2=np.floor(np.random.uniform(0,255,[64, 64]))
print('\n segnale I1 \n',I1)
print('\n segnale I2 \n',I2)

prod_scal=prodscal(I1,I2)
print('\n il risultato del prodotto scalare dei due segnali è \n',prod_scal)

prod_scal2=prodscal(I2,I1)  # stavolta cambio l'ordine dei segnali
if prod_scal==prod_scal2:
    print('\n la proprietà commutativa è verificata \n')
else:
    print('\n la proprietà commutativa non è verificata \n')

b=random.randint(-1000, 1000)  # scalare casuale
molt1=prodscal(b*I1,I2)   # prodotto fra scalare e I1
molt2=prodscal(I1,b*I2)
if molt1==molt2:
    print('\n la proprietà di omogeneità è verificata \n' )
else :  
    print('\n la proprietà di omogeneità è verificata \n' )

I3=np.floor(np.random.uniform(0,255,[64, 64]))  # terzo segnale
if prodscal(I1+I2,I3)==prodscal(I1,I3)+prodscal(I2,I3):
    if prodscal(I1,I2+I3)==prodscal(I1,I2)+prodscal(I1,I3):
        print('\n la proprietà distributiva è verificata \n')
    else:
        print('\n la proprietà distributiva non è verificata \n')
            



        


