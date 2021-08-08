from stuff import *
import numpy as np
from math import exp,log,sqrt
from scipy.optimize import minimize
from scipy.special import gammaln

mindate='2021-04-15'

# Alpha, Delta counts by day
A=np.zeros(1000,dtype=int)
D=np.zeros(1000,dtype=int)
with open('alphadelta','r') as fp:
  maxd=0
  for x in fp:
    y=x.strip().split()
    if y[0]>=mindate:
      d=datetoday(y[0])-datetoday(mindate)
      A[d]=int(y[2])
      D[d]=int(y[3])
      maxd=max(maxd,d)
ndays=maxd+1

def LL(xx):
  t0,lam=xx
  LL=0
  for day in range(ndays):
    e=exp(lam*(day-t0))
    a,d=A[day],D[day]
    if a==0 and d==0: continue
    LL+=a*log(1/(1+e))+d*log(e/(1+e))+(a+d)*(log(a+d)-1)-gammaln(a+1)-gammaln(d+1)
  return LL

condition=1e3
def NLL(xx): return -LL(xx)/condition

def Fisher(xx,eps=1e-4):
  t0,lam=xx
  fi=(-LL([t0,lam-eps])+2*LL([t0,lam])-LL([t0,lam+eps]))/eps**2
  zconf=1.96
  return zconf/sqrt(fi)

d0=datetoday('2021-05-15')-datetoday(mindate)
res=minimize(NLL,[0,0.1],bounds=[(d0-10,d0+10), (-0.2,0.2)], method="SLSQP")
if not res.success: raise RuntimeError(res.message)
t0,lam=res.x
dlam=Fisher(res.x)
print("%6.2f    %5.3f (%5.3f - %5.3f)     %8.3f"%(t0,lam,lam-dlam,lam+dlam,LL(res.x)))
