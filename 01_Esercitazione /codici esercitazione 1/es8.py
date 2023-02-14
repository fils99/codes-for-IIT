import numpy as np
import math
import matplotlib.pylab as plt

x=np.linspace(-5,5,100)
x0=float(input('immettere valore di x0:    '))
A=float(input('immettere valore di A:    '))
w=float(input('immettere valore di w:    '))
y=np.zeros(100)
for i in range(0,100):
    if x[i]>x0:
        if x[i]<x0+w:
            y[i]=A
    else:
        y[1]=0
fig, ax=plt.subplots()
plt.plot(x, y,'r--')
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.grid(True)
# plt.show() tanto non cambia nulla
               
    


