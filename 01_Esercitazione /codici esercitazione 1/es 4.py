import numpy as np
import math
import matplotlib.pylab as plt
from numpy import linalg as la
A=np.array([[1,2,3],[4,5,6],[7,8,9]])
B=np.diag(np.ones(3))
b=np.array([[3],[5],[7]])
som=A+A
sottr=A-B
prod=A*B
quadr=A**2
tens1=np.kron(A,B)
tens2=np.kron(B,A)
prod2=np.dot(A,b)
prod3=np.dot(np.transpose(b),A)