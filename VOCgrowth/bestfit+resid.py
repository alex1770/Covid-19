# py bestfit2+resid.py < file
# Format: Date, #Oldvariant, #Newvariant
# log(#New/#Old) ~= a + b*day
# Square errors are weighted by 1/(1/#Old+1/#New).

from stuff import *
import numpy as np
import sys
from math import sqrt

# Whether to interpret the first number as the total number of cases, or just the number of old variant cases
totalmode=(len(sys.argv)>1)

A=[];D=[];dt=[]
for x in sys.stdin:
  if x[0]=='#' or len(x.strip())==0: continue
  y=x.split()
  a,d=float(y[1]),float(y[2])
  if a>0 or d>0:
    dt.append(y[0])
    if totalmode: a-=d
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
mult=tres/n
print("Residual multiplier = %.3f"%mult)
# Imagining rescaling: A/=mult, D/=mult (so W/=mult, m/=mult, r/=mult) to rescale average residual to 1
# m/mult = -(Rescaled Hessian) of log likelihood = Observed Fisher information of rescaled problem
# Bottom right of inverse is 
gerr=sqrt(mult*np.linalg.solve(m,[0,1])[1])*1.96
print("Continuous growth rate = %.4f/day (%.4f/day - %.4f/day)"%(c[1],c[1]-gerr,c[1]+gerr))
d=-c[0]/c[1]
d1=int(round(-c[0]/c[1]))
print("Crossover on",daytodate(datetoday(dt[0])+d1)," %+.0f hours"%((d-d1)*24))

#A/=mult
#D/=mult
#W/=mult
#m/=mult
#r/=mult
#print(sqrt(np.linalg.solve(m,[0,1])[1])*1.96)
