#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 00:07:43 2022

@author: filippopasserini
"""

import numpy as np

""" esercizio 1 """
def DFTinv(matrice_coef,vettore_trasf,lung):
    DFTinv=np.dot(matrice_coef,vettore_trasf)
    for n in range(0,lung):
        DFTinv[n]=float(DFTinv[n])
    return DFTinv

""" esercizio 2 """

def prodscal(segnale1,segnale2):
    prodscal=0
    for n in range (0,64):
        for m in range (0,64):
            prodscal+=segnale1[n,m]*segnale2[n,m]
    return prodscal

    
