from stuff import *
import numpy as np
from math import exp,log,sqrt
import scipy
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
      A[d]=int(y[2])
      D[d]=int(y[3])
      maxd=max(maxd,d)
ndays=maxd+1

# LLt = log( P(NB(mu,p)=a) * P(NB(e*mu^rho,p)=d) )
def LLt(mu,a,d,e,rho,p):
  muA=mu*(1-p)/p
  muD=e*mu**rho*(1-p)/p
  LLa=gammaln(a+muA)-gammaln(muA)-gammaln(a+1)+muA*log(1-p)+a*log(p)
  LLd=gammaln(d+muD)-gammaln(muD)-gammaln(d+1)+muD*log(1-p)+d*log(p)
  return LLa+LLd

# LLt1 = d/dmu ( log( P(NB(mu,p)=a) * P(NB(e*mu^rho,p)=d) ) )
def LLt1(mu,a,d,e,rho,p):
  muA=mu*(1-p)/p;muA1=(1-p)/p
  muD=e*mu**rho*(1-p)/p;muD1=e*rho*mu**(rho-1)*(1-p)/p
  LLa1=(digamma(a+muA)-digamma(muA)+log(1-p))*muA1
  LLd1=(digamma(d+muD)-digamma(muD)+log(1-p))*muD1
  return LLa1+LLd1

def NLLt(xx,a,d,e,rho,p):
  return -LLt(xx[0],a,d,e,rho,p)

def NLLt1(xx,a,d,e,rho,p):
  return -LLt1(xx[0],a,d,e,rho,p)

# rho = T_a/T_d
# q = p/(1-p) = mu/r
def LL(xx):
  t0,lam,rho,q=xx
  p=q/(1+q)
  LL=0
  for day in range(ndays):
    e=exp(lam*(day-t0))
    a,d=A[day],D[day]
    if a==0 and d==0: continue

    mu0=1e-5
    mu1=max(-a*q/log(1-p),(-d*q/log(1-p)/e)**(1/rho))

    # (Annoyingly LLt isn't log-concave in general)
    assert LLt1(mu0,a,d,e,rho,p)>0 and LLt1(mu1,a,d,e,rho,p)<0
    while mu1-mu0>1e-12*(mu0+mu1):
      mu=(mu0+mu1)/2
      if LLt1(mu,a,d,e,rho,p)>0: mu0=mu
      else: mu1=mu
    mu=(mu0+mu1)/2
    
    LL+=LLt(mu,a,d,e,rho,p)
  return LL

condition=1e3
def NLL(xx): return -LL(xx)/condition

def Fisher(xx,eps=1e-4):
  t0,lam,rho,q=xx
  fi=(-LL([t0,lam-eps,rho,q])+2*LL([t0,lam,rho,q])-LL([t0,lam+eps,rho,q]))/eps**2
  zconf=1.96
  return zconf/sqrt(fi)


res=minimize(NLL,[50,0.1,1,1],bounds=[(40,60), (0,0.2), (0.5,2), (1e-4,100)], method="SLSQP")
if not res.success: raise RuntimeError(res.message)
t0,lam,rho,q=res.x
dlam=Fisher(res.x)
print("%6.2f    %5.3f (%5.3f - %5.3f)     %8.3f"%(t0,lam,lam-dlam,lam+dlam,LL(res.x)),rho,q)
