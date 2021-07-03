from stuff import *
from scipy.optimize import minimize
from scipy.special import gammaln
import numpy as np
import sys
from math import log,exp,sqrt,sin,pi
from random import random

#np.set_printoptions(precision=3,linewidth=120)


### Model ###
#
# Given
#
# n = number of steps (e.g., years of age)
# constr = {(a,b): c, ...}
# f()  (= identity or exp())
# 
# seek x[0],...,x[n-1] such that
#
# for all constr[(a,b)]=c, sum_{a<=i<b} f(x[i]) ~= c
# x[] is smooth
#
# Actually parametrise x[0..n-1] by X[0..N-1], via
#
# bmL = L = numsteps + 2*bmsig
# x_i = A+B*sqrt(L)*(i/L*X_0 + sqrt(2)/pi*sum_n e^{-(n*bmsig/L)^2/2}sin(n*pi*i/L)*X_n/n)
# where X_n ~ N(0,1),  n=0,...,N-1; N=ceil(3*L/bmsig) (say)
# N=bmN is just an approximation to infinity - no particular significance
#
### End Model ###


def expand(xx):
  w=xx[1]*sqrt(bmN)
  yy=[]
  for i in range(nsteps):
    t=i/bmL
    x=xx[0]+w*(t*xx[2]+np.dot(bmsin2[i],xx[2:]))
    yy.append(f(x))
  return np.array(yy)

# Return negative log likelihood (negative because scipy can only minimise, not maximise)
# If const is true then add in all the constant terms (that don't affect the optimisation)
def NLL(xx,const=False):
  tot=0
  
  yy=expand(xx)
  for (a,b) in constr:
    mu=yy[a:b].sum()
    n=constr[(a,b)]
    if mu==0: print(xx);print(yy)
    tot+=max((-mu+n*log(nif*mu))*nif,-10000)
    if const: tot+=log(nif)-gammaln(nif*n+1)# Approx normalisation
  
  for i in range(2+bmN):
    if i==0: var=100
    elif i==1: var=100
    else: var=1
    tot+=-xx[i]**2/(2*var)
    if const: tot-=log(2*pi*var)/2

  return -tot

#def f(x):
#  if x>0: return x+log(1+exp(-x))+1e-3
#  else: return log(1+exp(x))+1e-3
  
def f(x):
  if x<100: return exp(x)
  else: return 1e30
  
def smooth(nsteps0,constr0,bmsig=13,nif0=.1):
  global nsteps,constr,bmL,bmN,bmsin,bmweight,bmsin2,nif
  nsteps=nsteps0
  constr=constr0
  nif=nif0
  bmL=nsteps+int(bmsig*2+0.999)# Add on bmsig*2 to eliminate periodicity effects
  bmN=int(2.5*bmL/bmsig+1)
  bmsin=[sin(r*pi/bmL) for r in range(2*bmL)]
  bmweight=[0]+[sqrt(2)/pi*exp(-(n*bmsig/bmL)**2/2)/n for n in range(1,bmN)]
  bmsin2=[np.array([bmweight[n]*bmsin[(i*n)%(2*bmL)] for n in range(bmN)]) for i in range(nsteps)]
  
  # 2+bmN parameters to be optimised: A, B, X_0, ..., X_{bmN-1}
  
  xx=[0,0]+[(random()-.5)*.1 for i in range(bmN)]
  best=(-1e30,)
  for it in range(3):
    res=minimize(NLL,xx,method="SLSQP",bounds=[(-20,20),(0,10)]+[(-3,3)]*bmN,options={"maxiter":1000,"eps":1e-4})
    LL=-NLL(res.x,const=True)
    if LL>best[0]: best=(LL,res.x)
    xx=res.x+[(random()-.5)*.1 for i in range(bmN+2)]
    xx[1]=abs(xx[1])
  #if not res.success: raise RuntimeError(res.message)
  if not res.success: print("Warning:",res.message,file=sys.stderr)
  return expand(best[1]),best[0]
    
if __name__=='__main__':
  testnsteps=95
  testconstr={
    (0,5): 433,
    (5,10): 1321,
    (10,15): 2416,
    (15,20): 2887,
    (20,25): 4193,
    (25,30): 2183,
    (30,35): 1538,
    (35,40): 1331,
    (40,45): 1248,
    (45,50): 1025,
    (50,55): 810,
    (55,60): 508,
    (60,65): 305,
    (65,70): 195,
    (70,75): 127,
    (75,80): 76,
    (80,85): 54,
    (85,90): 36,
    (90,95): 10,
  }
  testconstr={(0, 5): 0, (5, 10): 0, (10, 15): 0, (15, 20): 0, (20, 25): 0, (25, 30): 1, (30, 35): 2, (35, 40): 2, (40, 45): 2, (45, 50): 2, (50, 55): 2, (55, 60): 4, (60, 65): 10, (65, 70): 5, (70, 75): 6, (75, 80): 9, (80, 85): 7, (85, 90): 18, (90, 150): 15}

  yy,LL=smooth(testnsteps,testconstr)

  with open('tempgr','w') as fp:
    for (i,y) in enumerate(yy): print(i,y,file=fp)
  print(LL)
  print(', '.join("%.2f"%y for y in yy))
  for (a,b) in constr:
    mu=yy[a:b].sum()
    n=testconstr[(a,b)]
    print("(%d,%d): %d   %g"%(a,b,n,mu))
  print()
