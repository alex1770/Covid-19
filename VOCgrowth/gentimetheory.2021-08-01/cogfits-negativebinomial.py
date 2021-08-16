from stuff import *
import numpy as np
from math import exp,log,sqrt
from scipy.optimize import minimize
from scipy.special import gammaln,digamma

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
      A[d]=int(y[1])
      D[d]=int(y[2])
      maxd=max(maxd,d)
ndays=maxd+1

def mumax(r,a,d,e,p):
  return gammaln(a+r)+gammaln(d+e*r)-gammaln(r)-gammaln(e*r)+r*(1+e)*log(1-p)

def mumaxd(r,a,d,e,p):
  return digamma(a+r)+e*digamma(d+e*r)-digamma(r)-e*digamma(e*r)+(1+e)*log(1-p)

# q = p/(1-p) = mu/r
def LL(xx):
  t0,lam,q=xx
  p=q/(1+q)
  LL=0
  for day in range(ndays):
    e=exp(lam*(day-t0))
    a,d=A[day],D[day]
    if a==0 and d==0: continue
    # Need to set r to maximise this:
    # Gamma(a+r)Gamma(d+e*r)/(Gamma(r)*Gamma(e*r))*(1-p)**(r*(1+e))
    r0=-1/((1+e)*log(1-p))
    r1=-1/((1+e)*log(1-p))*(a+d)
    assert mumaxd(r0,a,d,e,p)>=0 and mumaxd(r1,a,d,e,p)<=0
    while r1-r0>1e-12*(r0+r1):
      r=(r0+r1)/2
      if mumaxd(r,a,d,e,p)>=0: r0=r
      else: r1=r
    r=(r0+r1)/2
    LL+=gammaln(a+r)   - gammaln(a+1) - gammaln(r)     + r*log(1-p)+a*log(p)
    LL+=gammaln(d+e*r) - gammaln(d+1) - gammaln(e*r) + e*r*log(1-p)+d*log(p)
  return LL

condition=1e3
def NLL(xx): return -LL(xx)/condition

def Fisher(xx,eps=1e-4):
  t0,lam,q=xx
  fi=(-LL([t0,lam-eps,q])+2*LL([t0,lam,q])-LL([t0,lam+eps,q]))/eps**2
  zconf=1.96
  return zconf/sqrt(fi)

d0=datetoday('2021-05-15')-datetoday(mindate)
res=minimize(NLL,[0,0.1,1],bounds=[(d0-10,d0+10), (-0.2,0.2), (1e-4,100)], method="SLSQP")
if not res.success: raise RuntimeError(res.message)
t0,lam,q=res.x
dlam=Fisher(res.x)
print("%6.2f    %5.3f (%5.3f - %5.3f)     %8.3f"%(t0,lam,lam-dlam,lam+dlam,LL(res.x)),q)
