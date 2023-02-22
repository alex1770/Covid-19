# Beta binomial regression
# Arose from thinking about analogy Poisson:binomial = negative binomial:betabinomial
# BB(n,alpha,beta) in WP notation
# Var / (what would be var of binomial) = (alpha+beta+n)/(alpha+beta+1)
# Analogous to binomial regression, we're compressing variant numbers v0(t) given the sum v0(t)+v1(t)
# So the parameters of the BB at each timestep are allowed to depend on v0(t)+v1(t), but not v0(t) separately.
# Let's say there is global overdispersion parameter, m.
# And at each timestep we have the odds ratio, rho, given by exp(a+b*t) or some such.
# Then we choose beta/alpha=rho and (alpha+beta+n)/(alpha+beta+1)=m and (perforce) n=v0(t)+v1(t).
# Hmm, it will go mental if it wants to make m<=1.
# alpha+beta=(n-m)/(m-1)
# alpha=(alpha+beta)/(1+rho)
# beta=(alpha+beta)*rho/(1+rho)

import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-c', '--col0',        type=int,default=1,    help="Column number of variant 0 (includes date column, count from 0)")
parser.add_argument('-d', '--col1',        type=int,default=2,    help="Column number of variant 1 (includes date column, count from 0)")
parser.add_argument('-f', '--mindate',     default="2019-01-01",  help="Min sample date of sequence")
parser.add_argument('-t', '--maxdate',     default="9999-12-31",  help="Max sample date of sequence")
parser.add_argument('-p', '--prlevel',     type=int,default=1,    help="Print level")
args=parser.parse_args()

from stuff import *
import numpy as np
from math import exp,log,sqrt
from scipy.optimize import minimize
from scipy.special import gammaln,betaln
from scipy.stats import norm, multivariate_normal as mvn

conf=0.95
zconf=norm.ppf((1+conf)/2)
maxmult=20
prlevel=args.prlevel
print("Using date range",args.mindate,"-",args.maxdate)

# Variant0, Variant1 counts by day
N0=[];N1=[];DT=[]
if 0:
  #fp=open('BA5vsCent','r')
  #fp=open('UK_BA.5*_BA.2.12.1_BA.4*','r')
  #fp=open('Japan','r')
  fp=open('UK_BA.5*_BA.2.75')
  #fp=open("UK_BA.2_BA.2.23")
else:
  fp=sys.stdin

for x in fp:
  if x[0]=='#': continue
  y=x.strip().split()
  if y[0]>=args.mindate and y[0]<=args.maxdate:
    d=datetoday(y[0])
    try:
      v0=float(y[args.col0])
    except ValueError:
      v0=0
    try:
      v1=float(y[args.col1])
    except ValueError:
      v1=0
    if v0+v1>=1: N0.append(v0);N1.append(v1);DT.append(d)

minday=min(DT)
maxday=max(DT)+1
ndays=maxday-minday
V0=np.zeros(ndays)
V1=np.zeros(ndays)
for (v0,v1,dt) in zip(N0,N1,DT):
  V0[dt-minday]=v0
  V1[dt-minday]=v1

def LL(xx,pr=0):
  a0,lam,mult=xx
  LL=0
  for day in range(ndays):
    n0,n1=V0[day],V1[day]
    rho=exp(a0+lam*(day-day0))# odds ratio
    n=n0+n1
    if n>mult:
      s=(n-mult)/(mult-1)
    else:
      s=1e-3
    al=s/(1+rho)
    be=s*rho/(1+rho)
    # Now BB(n,al,be) has sample space {0,1,...,n}, mean=n/(1+rho) and overdispersion=mult
    DLL=gammaln(n0+n1+1)-gammaln(n0+1)-gammaln(n1+1)+betaln(n0+al,n1+be)-betaln(al,be)
    if pr:
      if n0>0 and n1>0: print(Date(day),"%6d %6d %7.3f %7.3f %10.2f %10.2f %7.3f"%(n0,n1,log(rho),log(n1/n0),al,be,DLL))
      else: print(Date(day),"%6d %6d %7.3f       - %10.2f %10.2f %7.3f"%(n0,n1,log(rho),al,be,DLL))
    LL+=DLL
  return LL

condition=1e3
def NLL(xx): return -LL(xx)/condition

def Hessian(xx,eps):
  N=len(xx)
  H=np.zeros([N,N])
  for i in range(N-1):
    for j in range(i+1,N):
      v=0
      for (s1,s2) in [(-1,-1),(-1,1),(1,-1),(1,1)]:
        x=np.copy(xx)
        x[i]+=s1*eps[i]
        x[j]+=s2*eps[j]
        v+=s1*s2*LL(x)
      e=v/(4*eps[i]*eps[j])
      H[i,j]=e
      H[j,i]=e
  for i in range(N):
    x=np.copy(xx)
    v=0
    for s in [-1,0,1]:
      x=np.copy(xx)
      x[i]+=s*eps[i]
      v+=(s*s*3-2)*LL(x)
    H[i,i]=v/eps[i]**2
  return H

def getCI(C):
  N=len(xx)
  err=[]
  for i in range(N):
    if C[i,i]>=0: err.append(zconf*sqrt(C[i,i]))
    else: err.append(None)
  return err

smooth=1
# Do simple regression to get decent initial values for NB regression
# V0s, V1s = smoothed V0, V1
V0s=np.zeros(len(V0))
V1s=np.zeros(len(V1))
for i in range(ndays):
  i0=max(i-smooth,0)
  i1=min(i+smooth+1,ndays)
  V0s[i]=V0[i0:i1].sum()/(i1-i0)
  V1s[i]=V1[i0:i1].sum()/(i1-i0)
V0sp=V0s+1e-30
V1sp=V1s+1e-30
W=1/(1/V0sp+1/V1sp)
day0=ndays//2
X=np.arange(ndays)-day0
Y=np.log(V1sp/V0sp)
m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
r=np.array([sum(W*Y),sum(W*X*Y)])
c=np.linalg.solve(m,r)
C=np.linalg.pinv(m)
dlam=zconf*sqrt(C[1,1])
if prlevel>=2: print("Simple regression growth: %.4f (%.4f - %.4f)  (but CI may be a bit off due to smoothing)"%(c[1],c[1]-dlam,c[1]+dlam))
rho=np.exp(c[0]+c[1]*X)
T=V0+V1
# Take off one degree of freedom because rho is tuned to V0s, V1s using two degrees of freedom, but
# we're evaluating residuals using the original V0, V1. (Semi-guess, with some empirical backup.)
mult0=((V1-V0*rho)**2/rho).sum()/(T.sum()-max(T))
if prlevel>=2: print("Variance overdispersion as estimated from simple residuals (though caution because smoothing): %.3f"%mult0)
# dayoffset=day0-c[0]/c[1]

desc=["Intercept","Growth of V1 rel V0","Overdispersion multiplier"]
ml=max(len(x) for x in desc)
bounds=[(c[0]-3,c[0]+3), (c[1]-0.2,c[1]+0.2), (1.01,maxmult)]
res=minimize(NLL,[c[0],c[1],min(2,maxmult)],bounds=bounds, method="SLSQP", options={'ftol':1e-20, 'maxiter':10000})
if not res.success: raise RuntimeError(res.message)
print("Log likelihood: %.3f"%(LL(res.x)))
a0,lam,mult=xx=res.x
eps=(5e-3,1e-4,1e-3)
H=Hessian(xx,eps)
negdef=np.all(np.linalg.eigvals(H)<0)
C=-np.linalg.inv(H)
da0,dlam,dmult=getCI(C)
for i in range(len(xx)):
  if xx[i]<bounds[i][0]+1e-3 or xx[i]>bounds[i][1]-1e-3: print("Warning:",desc[i],"hit bound")

print("Growth of V1 rel V0: %.4f"%lam,end="")
if dlam==None: print(" (couldn't evaluate CI)")
else: print(" (%.4f - %.4f)"%(lam-dlam,lam+dlam))

print("Crossover date: %s"%(Date(minday+int(round(day0-a0/lam)))),end="")
if not negdef:
  print(" (CI n/a as Hessian is not negative definite)")
elif lam-dlam<1e-3 and lam+dlam>-1e-3:
  print(" (CI n/a as relative growth can get too close to 0)")
else:
  # Sampling is most convenient way of getting CrIs for crossover points
  numsamp=100000
  test=mvn.rvs(mean=xx,cov=C,size=numsamp)
  cr0=Date(minday+int(round(day0-np.quantile(test[:,0]/test[:,1],(1+conf)/2))))
  cr1=Date(minday+int(round(day0-np.quantile(test[:,0]/test[:,1],(1-conf)/2))))
  print(" (%s - %s)"%(cr0,cr1))

print("Variance overdispersion (proper estimate): %.3f"%mult,end="")
if dmult==None: print(" (couldn't evaluate CI)")
else: print(" (%.3f - %.3f)"%(mult-dmult,mult+dmult))

if prlevel>=3:
  print()
  print("                                                  =Smoothed log odds ratios=")
  print("Date       Num V0  Num V1       Pred     Actual   Pred_lor  Act_lor    Resid")
s0=s1=0
for day in range(ndays):
  a,d=V0[day],V1[day]
  if a==0 and d==0: continue
  rho=exp(a0+lam*(day-day0))
  #rho=exp(c[0]+c[1]*(day-day0))
  pred_pr=rho/(1+rho)
  actual_pr=d/(a+d)
  if prlevel>=3:
    print(Date(minday+day),"%6d  %6d  %9.5f  %9.5f"%(a,d,pred_pr,actual_pr),end="")
  a_s,d_s=V0s[day],V1s[day]
  if a_s>0 and d_s>0:
    pred_lor=log(rho)
    actual_lor=log(d_s/a_s)
    iv=a_s*d_s/(a_s+d_s)
    s0+=1
    s1+=(pred_lor-actual_lor)**2*iv
    if prlevel>=3: print("    %7.3f  %7.3f  %7.3f"%(pred_lor,actual_lor,(pred_lor-actual_lor)*sqrt(iv)))
  else:
    if prlevel>=3:
      print()

#print("Variance overdispersion as estimated from BB residuals (though caution because smoothing): %.3f"%(s1/s0))
