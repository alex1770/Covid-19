import numpy as np
from math import log,exp
from random import random

def getdiffmat(n,order):
  A=np.zeros([n,n])
  if order==1:
    for i in range(n-1):
      A[i,i]+=1
      A[i+1,i+1]+=1
      A[i,i+1]-=1
      A[i+1,i]-=1
  elif order==2:
    for i in range(n-2):
      A[i,i]+=1
      A[i,i+1]+=-2
      A[i,i+2]+=1
      A[i+1,i]+=-2
      A[i+1,i+1]+=4
      A[i+1,i+2]+=-2
      A[i+2,i]+=1
      A[i+2,i+1]+=-2
      A[i+2,i+2]+=1
  else: raise RuntimeError("Unrecognised order %d"%order)
  return A

order=2
n=100

A=getdiffmat(n,2)
for i in range(9):
  for j in range(i+1,10):
    B=np.delete(np.delete(A,[i,j],0),[i,j],1)
    print(i,j,(j-i)**2,int(round(np.linalg.det(B))))

def detpert(eps,order):
  if order==1:
    return eps.sum()
  elif order==2:
    n=len(eps)
    ind=np.arange(n)
    s0=eps.sum()
    s1=(eps*ind).sum()
    s2=(eps*ind*ind).sum()
    return s0*s2-s1*s1
  else:
    assert 0
    
for iv in [1,10,100,1000,10000,1e5,1e6,1e7,1e8,1e9,1e10]:
  #  Q(x[], al[], iv) = sum_i { eps[i].x[i]^2 + iv*sum_i (x[i+1]-x[i])^2 }
  #  lim_{eps[]->0} (sum(eps[]))^{1/2}*\int_{x[]} exp(-(1/2)Q(x[], eps[], iv)) = iv^{-(1/2)(n-1)}
  for eps in [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]:
    A=np.zeros([n,n])
    epsarray=np.array([eps*iv*random() for i in range(n)])
    for i in range(n):
      A[i,i]+=epsarray[i]
    A+=iv*getdiffmat(n,order)
    M1=-(1/2)*np.linalg.slogdet(A)[1]+(1/2)*log(detpert(epsarray,order))
    N1=-(1/2)*(n-order)*log(iv)
    print(eps,M1,N1,M1-N1)
  print()
  
