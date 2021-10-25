from stuff import *
import sys
import numpy as np
from math import exp,log,sqrt
from scipy.optimize import minimize
from scipy.special import gammaln

mindate='2021-04-20'
maxdate='2021-05-31'
if len(sys.argv)>1: mindate=sys.argv[1]
if len(sys.argv)>2: maxdate=sys.argv[2]
print("Using date range",mindate,"-",maxdate)

# Alpha, Delta counts by day
A=np.zeros(1000,dtype=int)
D=np.zeros(1000,dtype=int)
#with open('alphadelta','r') as fp:
fp=sys.stdin
if 1:
  maxd=0
  for x in fp:
    y=x.strip().split()
    if y[0]>=mindate and y[0]<=maxdate:
      d=datetoday(y[0])-datetoday(mindate)
      A[d]=int(y[1])
      D[d]=int(y[2])
      maxd=max(maxd,d)
ndays=maxd+1

def LL(xx):
  t0,lam=xx
  LL=0
  for day in range(ndays):
    e=exp(lam*(day-t0))
    a,d=A[day],D[day]
    if a==0 and d==0: continue
    # max_{r,s | s=r*e} log(P(Po(r)=a and Po(s)=d))
    LL+=a*log(1/(1+e))+d*log(e/(1+e))+(a+d)*(log(a+d)-1)-gammaln(a+1)-gammaln(d+1)
  return LL

condition=1e3
def NLL(xx): return -LL(xx)/condition

def Fisher(xx,eps=1e-4):
  t0,lam=xx
  fi=(-LL([t0,lam-eps])+2*LL([t0,lam])-LL([t0,lam+eps]))/eps**2
  zconf=1.96
  return zconf/sqrt(fi)

minday=datetoday(mindate)
d0=datetoday('2021-05-15')-minday
res=minimize(NLL,[d0,0.1],bounds=[(d0-20,d0+20), (-0.2,0.2)], method="SLSQP")
if not res.success: raise RuntimeError(res.message)
t0,lam=res.x
dlam=Fisher(res.x)
print("%s  %6.2f    %6.4f (%6.4f - %6.4f) %7.5f    %8.3f"%(daytodate(minday+int(round(t0))),t0,lam,lam-dlam,lam+dlam,dlam,LL(res.x)))
