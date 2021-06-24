# py bestfit.py < alphadelta
# Fits cols 3, 4 (Alpha, Delta) using
# log(Delta/Alpha) ~= a + b*day
# Square errors are weighted by 1/(1/a+1/b).

from stuff import *
import numpy as np
import sys

A=[];D=[];dt=[]
for x in sys.stdin:
  y=x.split()
  dt.append(y[0])
  A.append(float(y[2]))
  D.append(float(y[3]))

n=len(A)
A=np.array(A)
D=np.array(D)
W=1/(1/A+1/D)
dt0=datetoday(dt[0])
X=[datetoday(d)-dt0 for d in dt]
Y=np.log(D/A)

m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
r=np.array([sum(W*Y),sum(W*X*Y)])
c=np.linalg.solve(m,r)
print("Continuous growth rate = %.3f/day"%(c[1]))
print("Crossover on",daytodate(datetoday(dt[0])+int(-c[0]/c[1]+.5)))
