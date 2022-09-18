from stuff import *
import numpy as np
from math import exp,log,sqrt
from scipy.optimize import minimize
from scipy.special import gammaln,digamma

mindate='0000-00-00'
maxdate='9999-99-99'
prlevel=1
if len(sys.argv)>1: mindate=sys.argv[1]
if len(sys.argv)>2: maxdate=sys.argv[2]
if len(sys.argv)>3: prlevel=int(sys.argv[3])
print("Using date range",mindate,"-",maxdate)

# Variant0, Variant1 counts by day
V0=[];V1=[];DT=[]
fp=sys.stdin
#fp=open('alphadelta','r')
if 1:
  for x in fp:
    y=x.strip().split()
    if y[0]>=mindate and y[0]<=maxdate:
      d=datetoday(y[0])
      v0=int(y[1])
      v1=int(y[2])
      if v0>0 and v1>0: V0.append(v0);V1.append(v1);DT.append(d)
ndays=len(V0)

def mumax(r,a,d,e,p):
  return gammaln(a+r)+gammaln(d+e*r)-gammaln(r)-gammaln(e*r)+r*(1+e)*log(1-p)

def mumaxd(r,a,d,e,p):
  return digamma(a+r)+e*digamma(d+e*r)-digamma(r)-e*digamma(e*r)+(1+e)*log(1-p)

# q = p/(1-p) = mu/r
def LL(xx):
  t0,lam,q=xx
  p=q/(1+q)
  LL=0
  for i in range(ndays):
    a,d,day=V0[i],V1[i],DT[i]
    e=exp(lam*(day-dayoffset-t0))
    # Setting r to maximise P( NB(r,p)=a and NB(er,p)=d ), which is a constant times this:
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

def Fisher(xx):
  zconf=1.96
  t0,lam,q=xx
  eps=(5e-3,1e-4,1e-3)
  H=Hessian(xx,eps)
  HI=np.linalg.inv(H)
  N=len(xx)
  err=[zconf*sqrt(-HI[i,i]) for i in range(N)]
  return err

# Do simple regression to get good initial values for NB regression
A=np.array(V0)
D=np.array(V1)
W=1/(1/A+1/D)
day0=int(round(sum(DT)/len(DT)))
X=np.array(DT)-day0
Y=np.log((D+1e-20)/(A+1e-20))
m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
r=np.array([sum(W*Y),sum(W*X*Y)])
c=np.linalg.solve(m,r)
mi=np.linalg.pinv(m)
dayoffset=day0-c[0]/c[1]
# Could estimate overdispersion parameter, but it's not likely to be that close to what NB finds because NB
# is quite insensitive to its exact value. Which actually suggests there is not much point in using NB instead
# of Poisson. Ho hum.
# R=c[0]+c[1]*X-Y
# q0=(R*R*W).sum()/len(R)-1

res=minimize(NLL,[0,c[1],1],bounds=[(-20,20), (-0.2,0.2), (1e-4,100)], method="SLSQP", options={'ftol':1e-20, 'maxiter':10000})
if not res.success: raise RuntimeError(res.message)
t0,lam,q=res.x
dt0,dlam,dq=Fisher(res.x)
print("Log likelihood: %.3f"%(LL(res.x)))
print("Growth of V1 rel V0: %.4f (%.4f - %.4f)"%(lam,lam-dlam,lam+dlam))
print("Crossover date: %s"%(daytodate(int(round(dayoffset+t0)))))
print("Variance overdispersion: %.3f (%.3f - %.3f)"%(1+q,1+q-dq,1+q+dq))
if prlevel>=2:
  print()
  print("Num V0  Num V1     Pred   Actual    Resid")
s0=s1=0
for i in range(ndays):
  a,d,day=V0[i],V1[i],DT[i]
  if a==0 or d==0: continue
  pred=lam*(day-dayoffset-t0)
  actual=log(d/a)
  v=1/a+1/d
  s0+=1
  s1+=(pred-actual)**2/v
  if prlevel>=2:
    print("%6d  %6d  %7.3f  %7.3f  %7.3f"%(a,d,pred,actual,(pred-actual)/sqrt(v*(1+q))))
print("Variance overdispersion as estimated from residuals: %.3f"%(s1/s0))
