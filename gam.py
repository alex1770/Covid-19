from stuff import *
from scipy.optimize import minimize
from scipy.special import gammaln
import numpy as np
from math import log,exp,sqrt,sin,pi

np.set_printoptions(precision=3,linewidth=120)


nsteps=95
constr={
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
bmsig=13
nif=1

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


bmL=nsteps+int(bmsig*2+0.999)# Add on bmsig*2 to eliminate periodicity effects
bmN=int(2.5*bmL/bmsig+1)
bmsin=[sin(r*pi/bmL) for r in range(2*bmL)]
bmweight=[0]+[sqrt(2)/pi*exp(-(n*bmsig/bmL)**2/2)/n for n in range(1,bmN)]
bmsin2=[np.array([bmweight[n]*bmsin[(i*n)%(2*bmL)] for n in range(bmN)]) for i in range(nsteps)]

# bmN parameters to be optimised: A, B, X_0, ..., X_{bmN-1}

def f(x):
  if x>0: return x+log(1+exp(-x))+1e-3
  else: return log(1+exp(x))+1e-3
  
def f(x):
  if x<100: return exp(x)
  else: return 1e30
  
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
def NLL(xx):
  tot=0
  
  yy=expand(xx)
  #print(xx)
  #print(yy)
  for (a,b) in constr:
    mu=yy[a:b].sum()
    n=constr[(a,b)]
    if mu==0: print(xx);print(yy)
    tot+=max((-mu+n*log(nif*mu))*nif,-10000)
    tot+=log(nif)-gammaln(nif*n+1)# Approx normalisation
  
  for i in range(bmN):
    tot+=-xx[i]**2/2
  tot-=bmN*log(2*pi)/2

  return -tot

xx=[1,1]+[0]*bmN
while 1:
  res=minimize(NLL,xx,method="SLSQP",bounds=[(0,20),(0,10)]+[(-3,3)]*bmN,options={"maxiter":1000,"eps":1e-4})
  xx=res.x
  print(res.message)
  yy=expand(xx)
  with open('tempgr','w') as fp:
    for (i,y) in enumerate(yy): print(i,y,file=fp)
  print(NLL(xx))
  print([int(y+.5) for y in yy])
  for (a,b) in constr:
    mu=yy[a:b].sum()
    n=constr[(a,b)]
    print("(%d,%d): %d   %g"%(a,b,n,mu))
  print()
  if res.success: break
