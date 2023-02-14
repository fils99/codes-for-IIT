import numpy as np
import math
import matplotlib.pylab as plt
a=np.arange(1,11,1)
b=np.array([[10],[9],[8],[7],[6],[5],[4],[3],[2],[1]])
c=np.arange(0,1.1,0.1)

d=np.array([-2, 4.5, 7])
f=np.array([0, -1, 3.33])
p=float(input('immettere un numero reale ≥ 1:    '))


def norma( vettore ):
    
      
    
    lung=len(vettore)
    
    norma_mom=0
    for j in range(0,lung):   
        norma_mom= norma_mom+(abs(vettore[j]))**p
    norma=float(norma_mom)**(1/p)
    return norma;

#norma dei vettori a, b, c
#norma_a=norma(a)
#norma_b=norma(b)
#norma_c=norma(c)

norma_d=norma(d)
norma_f=norma(f)

#norma di un vettore costruito a piacere
#lung_2=int(input('immettere lunghezza vettore (numero intero>0):   '))
#vettore_costruito=[]
#for i in range(0,lung_2):
#    el=input('immettere componente vettore (numero floating):   ')
#    vettore_costruito.append(float(el))
#norma_v=norma(vettore_costruito)

norma_somma=norma(d+f)
norma_sott=norma(d-f)
primo_membro=norma_somma**p+norma_sott**p
secondo_membro=p*(norma(d)**p+norma(f)**p)
