import numpy as np
import math
from numpy import linalg as la
t=np.arange(-1,1.01,0.01)

# scegliere una coppia di a e b

#a=np.sin(2*math.pi*t)
#b=t

#a=np.sin(6*math.pi*t)
#b=np.sin(4*math.pi*t)

# a=np.cos(3/2*math.pi*t) 
# b=t

a=np.exp(4*math.pi*1j*t)
b=np.exp(4*math.pi*1j*t)

b_con=np.conjugate(b)
import matplotlib.pylab as plt
plt.plot(t,a)
plt.plot(t,b)
prod=np.inner(a,b_con)