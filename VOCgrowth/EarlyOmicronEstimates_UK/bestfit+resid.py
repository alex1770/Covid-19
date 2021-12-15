# py bestfit2+resid.py < file
# Format: Date, #Oldvariant, #Newvariant
# log(#New/#Old) ~= a + b*day
# Square errors are weighted by 1/(1/#Old+1/#New).

from stuff import *
import numpy as np
import sys
from math import sqrt

A=[];D=[];dt=[]
for x in sys.stdin:
  y=x.split()
  a,d=float(y[1]),float(y[2])
  if a>0 or d>0:
    dt.append(y[0])
    A.append(a)
    D.append(d)

n=len(A)
A=np.array(A)
D=np.array(D)
W=A*D/(A+D)
dt0=datetoday(dt[0])
X=np.array([datetoday(d)-dt0 for d in dt])
Y=np.log((D+1e-30)/(A+1e-30))

m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
r=np.array([sum(W*Y),sum(W*X*Y)])
c=np.linalg.solve(m,r)
res=c[0]+c[1]*X-Y
tres=(W*res*res).sum()
print("Continuous growth rate = %.4f/day"%(c[1]))
d=-c[0]/c[1]
d1=int(round(-c[0]/c[1]))
print("Crossover on",daytodate(datetoday(dt[0])+d1)," %+.0f hours"%((d-d1)*24))
print("Residual multiplier = %.3f"%sqrt(tres/n))
