# https://jax.readthedocs.io/en/latest/notebooks/quickstart.html

from jax import jacfwd, jacrev
import jax.numpy as jnp
from jax import grad, jit, vmap
from random import random
from math import log,exp,sqrt
from scipy.stats import norm, multivariate_normal as mvn
from scipy.optimize import minimize
from jax.scipy.optimize import minimize as jminimize
from time import process_time as cpu
import numpy as np
from jax.config import config
config.update("jax_enable_x64", True)
jnp.set_printoptions(precision=6,suppress=True,linewidth=200)
conf=0.95

# Odds, o = exp(x0 + x1*mpg + x2*(cyl==6) + x3*(cyl==8))
# Prob(am=1) = o/(1+o)

data=[]
with open('cardata') as fp:
  for x in fp:
    y=x.strip()
    if y=='' or y[0]=='#': continue
    z=y.split()
    # mpg (float), cyl (4, 6 or 8), am (0 or 1)
    data.append([float(z[0]),int(z[1]),int(z[8])])

data_in=[];data_out=[]
for (mpg,cyl,am) in data:
  data_in.append([1,mpg,cyl==6,cyl==8])
  data_out.append(am)# Try float too
data_in=jnp.array(data_in)# Input: 1, mpg, cyl==6, cyl==8
data_out=jnp.array(data_out)# Output: am

def NLL0(xx):
  LL=0
  for (mpg,cyl,am) in data:
    lp=xx[0]+xx[1]*mpg+xx[2]*(cyl==6)+xx[3]*(cyl==8)
    o=exp(lp)
    p=(am*o+1-am)/(1+o)
    LL+=log(p)
  return -LL

data_in_np=np.array(data_in)
data_out_np=np.array(data_out)
def NLL1(xx):
  oo=np.exp(data_in_np@xx)
  return -(np.log((1-data_out_np*(1-oo))/(1+oo))).sum()

@jit
def NLL(xx):
  oo=jnp.exp(data_in@xx)
  return -(jnp.log((1-data_out*(1-oo))/(1+oo))).sum()

d1=jit(grad(NLL))
d2=jit(jacfwd(jacrev(NLL)))

# mean = -8.342991591730756, 0.3699446195026917, 0.7320799305429727, 0.701955898086171

bounds=[(-10,10), (0,1), (0,1), (0,1)]
t0=cpu()
res=minimize(NLL,[0.,0.5,0.5,0.5],jac=d1,bounds=bounds, method="SLSQP")
print(res.message)
print(cpu()-t0,"seconds")
print(res.fun)
print(res.x)
print()

# Seems not to work
# res=jminimize(NLL,jnp.array([0.,0.5,0.5,0.5]),method="BFGS")

if 1:
  t0=cpu()
  #xx=jnp.array([-8, 0.4, 0.7, 0.7])
  xx=res.x
  for it in range(20):
    ee=jnp.linalg.solve(d2(xx),d1(xx))
    xx-=ee
    #print(cpu()-t0,ee,xx)

xx=res.x
hes=d2(xx)
C=jnp.linalg.inv(hes)

# Lazy CI with MVN
# (Could actually get this exactly: lp will be a single normal, then use Phi as usual)
# Newdata prob = 0.3602787304187247
# CI 0.07257582754622581 - 0.8020835186992356
# CI 0.07259068495340207 - 0.8021377063775462
numsamp=1000000
print("Sampling")
samp=mvn.rvs(mean=xx,cov=C,size=numsamp)
print("Sampled")
newdata=[1, 21, 0, 0]
o=exp(xx@newdata)
print("Newdata prob =",o/(1+o))
lp=samp@newdata
lp.sort()
n0=int((1-conf)/2*numsamp)
n1=int((1+conf)/2*numsamp)
lp=lp[[n0,n1]]
oo=np.exp(lp)
pp=oo/(1+oo)
print("Normal CI from",numsamp,"MVN samples",pp[0],"-",pp[1])

# Doing the normal CI exactly
# Exact normal CI 0.07256838553436519 - 0.8021173166226204
mu=newdata@xx
v=newdata@np.array(C)@newdata
zconf=norm.ppf((1+conf)/2)
# X=N(mu,v); P(exp(X)/(1+exp(X)) < p_low) = (1-conf)/2
# l=mu-zconf*sqrt(v)
# exp(l)/(1+exp(l))
p0=1/(1+exp(-(mu-zconf*sqrt(v))))
p1=1/(1+exp(-(mu+zconf*sqrt(v))))
print("Exact normal CI",p0,"-",p1)
print()

# Now try MCMC over xx, with uniform prior
#        it       acc      time/it    CI
# 351137232  0.262010    220.245us    0.04863 - 0.80503
it=ac=0
al=2
L=-NLL(xx)
t0=cpu()
l=[]
nextprint=1000
while 1:
  yy=mvn.rvs(mean=xx,cov=al*C)
  L1=-NLL1(yy)
  if random()<exp(L1-L):
    xx=yy;L=L1
    ac+=1
  l.append(xx@newdata)
  it+=1
  if it==nextprint:
    print("%10d  %8.6f   %8.3fus"%(it,ac/it,(cpu()-t0)/it*1e6),end="")
    l.sort()
    n0=int((1-conf)/2*it)
    n1=int((1+conf)/2*it)
    lp=np.array([l[n0],l[n1]])
    oo=np.exp(lp)
    pp=oo/(1+oo)
    print("    %7.5f - %7.5f"%tuple(pp))
    nextprint=int(1.1*nextprint)
