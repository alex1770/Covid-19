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

def mumaxd(mu,a,d,e,rho):
  return mu+e*rho*mu**rho-(a+rho*d)

# rho = T_a/T_d
def LL(xx):
  t0,lam,rho=xx
  LL=0
  for day in range(ndays):
    e=exp(lam*(day-t0))
    a,d=A[day],D[day]
    if a==0 and d==0: continue
    # Need to set mu to maximise this: -mu-e*mu**rho+(a+d*rho)*log(mu)
    # which means solving this: mu+e*rho*mu**rho = a+rho*d
    mu0=0
    mu1=a+rho*d
    assert mumaxd(mu0,a,d,e,rho)<=0 and mumaxd(mu1,a,d,e,rho)>=0
    while mu1-mu0>1e-12*(mu0+mu1):
      mu=(mu0+mu1)/2
      if mumaxd(mu,a,d,e,rho)<=0: mu0=mu
      else: mu1=mu
    mu=(mu0+mu1)/2
    LL+=-mu+a*log(mu)-gammaln(a+1)
    LL+=-e*mu**rho+d*log(e*mu**rho)-gammaln(d+1)
  return LL

condition=1e3
def NLL(xx): return -LL(xx)/condition

def Fisher(xx,eps=1e-4):
  t0,lam,q=xx
  fi=(-LL([t0,lam-eps,q])+2*LL([t0,lam,q])-LL([t0,lam+eps,q]))/eps**2
  zconf=1.96
  return zconf/sqrt(fi)

res=minimize(NLL,[0,0.1,1],bounds=[(40,60), (-0.2,0.2), (0.1,10)], method="SLSQP")
if not res.success: raise RuntimeError(res.message)
t0,lam,rho=res.x
dlam=Fisher(res.x)
print("%6.2f    %5.3f (%5.3f - %5.3f)     %8.3f"%(t0,lam,lam-dlam,lam+dlam,LL(res.x)),rho)
