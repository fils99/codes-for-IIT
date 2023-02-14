import numpy as np
import math
from numpy import linalg as la
l1=int(input('immettere un numero intero ≥ 1:    '))
l2=int(input('immettere un numero intero ≥ 1:    '))
vettore_costruito1=[]
for i in range(0,l1):
    el=input('immettere componente vettore 1 (numero floating):   ')
    vettore_costruito1.append(complex(el))
vettore_costruito2=[]
for i in range(0,l2):
    el=input('immettere componente vettore 2 (numero floating):   ')
    vettore_costruito2.append(complex(el))
if l1==l2:
    def prodotto_scalare (vettore1,vettore2):
               prodotto_scalare=0
               for i in range(0,l1):
                   prodotto_scalare+=vettore1[i]*np.conj(vettore2[i])
               return prodotto_scalare;
    prodotto_scalare1_2=prodotto_scalare(vettore_costruito1,vettore_costruito2)  
else:
     print('errore - la lunghezza dei vettori non è compatibile ')
     print('non è possibile calcolare il prodotto scalare')                                    
# bo=np.inner(vettore_costruito1,np.conj(vettore_costruito2))  per verificare di aver fatto bene il prodotto                                                                                          