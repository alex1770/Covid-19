# https://jax.readthedocs.io/en/latest/notebooks/quickstart.html

import sys
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
jnp.set_printoptions(precision=6,suppress=True,linewidth=1000)
conf=0.95
ncat=3

# o_i = Odds(cyl==[4,6,8][i]) = exp(x_{i0}*1 + x_{i1}*mpg + x_{i2}*am)
# Prob(cyl==[4,6,8][i]) = o_i/(o_0+o_1+o_2)
# May gauge fix x_{00}=x_{01}=x_{02}=0

data=[]
with open('cardata') as fp:
  for x in fp:
    y=x.strip()
    if y=='' or y[0]=='#': continue
    z=y.split()
    # mpg (float), cyl (4, 6 or 8), am (0 or 1)
    data.append([float(z[0]),int(z[1]),int(z[8])])

nsmpg2=[]
with open('nsmpg2') as fp:
  for x in fp:
    nsmpg2.append([float(z) for z in x.strip().split()])
    
nsmpg1=[]
with open('nsmpg1') as fp:
  for x in fp:
    nsmpg1.append([float(z) for z in x.strip().split()])

# Modification to keep it bounded for df=0,1
data[4][0]=data[18][0]
nsmpg1[4]=nsmpg1[18]
nsmpg2[4]=nsmpg2[18]
# Extra modification required for df=2
data[2][0]=data[6][0]
nsmpg1[2]=nsmpg1[6]
nsmpg2[2]=nsmpg2[6]

# df=0 means use mpg directly
# df=1 means use ns(mpg,df=1)
# df=2 means use ns(mpg,df=2)
df=2
mult=1
if len(sys.argv)>1: mult=int(sys.argv[1])
print("Using multiplicity",mult)

data_in=[];data_out=[];data_out_cat=[]
for (i,(mpg,cyl,am)) in enumerate(data):
  if df==0:
    data_in.append([1,mpg,am])
  elif df==1:
    data_in.append([1,nsmpg1[i][0],am])
  elif df==2:
    data_in.append([1,nsmpg2[i][0],nsmpg2[i][1],am])
  else: assert 0
  data_out.append([float(cyl==4),float(cyl==6),float(cyl==8)])
  data_out_cat.append(int((cyl==6)+2*(cyl==8)))
data_in=jnp.array(data_in)# Input: 1, mpg or ns(mpg), am
data_out=jnp.array(data_out)# Output: Cat(cyl)
nvar=data_in.shape[1]

def getyop(xx):
  yy=np.zeros([ncat,nvar])
  yy[1:,:]=xx.reshape([ncat-1,nvar])
  oo=np.exp(data_in_np@yy.T)
  pp=((oo/oo.sum(axis=1)[:,None])*data_out_np).sum(axis=1)
  return yy,oo,pp

data_in_np=np.array(data_in)
data_out_np=np.array(data_out)
def NLL1(xx):
  yy=np.zeros([ncat,nvar])
  yy[1:,:]=xx.reshape([ncat-1,nvar])
  oo=np.exp(data_in_np@yy.T)
  pp=((oo/oo.sum(axis=1)[:,None])*data_out_np).sum(axis=1)
  return -(np.log(pp)).sum()*mult

@jit
def NLL(xx):
  yy=jnp.concatenate([jnp.zeros([1,nvar]),xx.reshape([ncat-1,nvar])])
  oo=jnp.exp(yy@data_in.T)
  pp=((oo/oo.sum(axis=0)[None,:])*data_out.T).sum(axis=0)
  return -(jnp.log(pp)).sum()*mult

d1=jit(grad(NLL))
d2=jit(jacfwd(jacrev(NLL)))

def boundscheck(xx,bb):
  for i in range(len(xx)):
    if xx[i]<bb[i][0]+1e-6 or xx[i]>bb[i][1]-1e-6: print("Variable",i,"=",xx[i],"hit bound",bb[i])

if df==0:
  bounds=[(-30,30), (-2,2), (-3,3)]*(ncat-1)
elif df==1:
  bounds=[(-30,30), (-100,100), (-3,3)]*(ncat-1)
elif df==2:
  bounds=[(-300,300), (-200,200), (-200,200), (-50,50)]*(ncat-1)
t0=cpu()
#res=minimize(NLL,[0.]*nvar*(ncat-1),bounds=bounds, method="SLSQP")
res=minimize(NLL,[0.]*nvar*(ncat-1),jac=d1,bounds=bounds, method="SLSQP")
boundscheck(res.x,bounds)
print(res.message)
print(cpu()-t0,"seconds")
xx=res.x
print(xx,NLL(xx))
print()

if 1:
  t0=cpu()
  xx=res.x
  for it in range(10):
    ee=jnp.linalg.solve(d2(xx),d1(xx))
    xx-=ee
    #print(cpu()-t0,ee,xx)

hes=d2(xx)
C=jnp.linalg.inv(hes)
print("NLL =",NLL(xx))
yy,oo,pp=getyop(xx)
print(yy)
print()

newdata=data_in[0]
o=np.exp(yy@newdata)
p=o/o.sum()
print("Newdata probs =",p)
numsamp=1000000

# CI with MVN (on xx, the log scale)
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26182 ( 0.03000 -  0.78178)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73749 ( 0.17680 -  0.96721)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00070 ( 0.00000 -  0.25222)
t0=cpu()
print("Sampling xx")
samp=mvn.rvs(mean=xx,cov=C,size=numsamp)
print("Sampled")
yysamp=jnp.concatenate([jnp.zeros([numsamp,1,nvar]),samp.reshape([numsamp,ncat-1,nvar])],axis=1)
oo=jnp.exp(yysamp@newdata)
pp=oo/oo.sum(axis=1)[:,None]
for i in range(ncat):
  print("Using %d samples of MVN: Prob(newdata = cat%d) = %8.5f (%8.5f - %8.5f)"%(numsamp,i,p[i],np.quantile(pp[:,i],(1-conf)/2),np.quantile(pp[:,i],(1+conf)/2)))
print("Time =",cpu()-t0)
print()

# Same, but do dot product with newdata first (should be identical results)
# df=2, mult=1
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26182 ( 0.03000 -  0.78180)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73749 ( 0.17678 -  0.96720)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00070 ( 0.00000 -  0.25237)
#
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26182 ( 0.14249 -  0.43042)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73749 ( 0.56784 -  0.85685)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00070 ( 0.00009 -  0.00508)
#
# df=2, mult=100
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26182 ( 0.21819 -  0.31067)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73749 ( 0.68848 -  0.78119)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00070 ( 0.00037 -  0.00131)
# Time = 26.5s

mu=xx.reshape([ncat-1,nvar])@newdata
D=(C.reshape([ncat-1,nvar,ncat-1,nvar])*newdata[None,:,None,None]*newdata[None,None,None,:]).sum(axis=[1,3])
t0=cpu()
print("Sampling xx@newdata")
samp=mvn.rvs(mean=mu,cov=D,size=numsamp)
print("Sampled")
oo=jnp.concatenate([jnp.ones([numsamp,1]),jnp.exp(samp)],axis=1)
pp=oo/oo.sum(axis=1)[:,None]
for i in range(ncat):
  print("Using %d samples of MVN: Prob(newdata = cat%d) = %8.5f (%8.5f - %8.5f)"%(numsamp,i,p[i],np.quantile(pp[:,i],(1-conf)/2),np.quantile(pp[:,i],(1+conf)/2)))
print("Time =",cpu()-t0)
print()

# CI with MVN on exp(xx@newdata) (kind of stupid, but including anyway because can do it semi-analytically)
# df=2, mult=1
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26182 (-0.29006 -  0.29257)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73749 (-2.79991 -  2.90877)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00070 (-2.03564 -  3.92320)
#
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26182 ( 0.15523 -  0.60348)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73749 ( 0.39408 -  0.84406)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00070 (-0.00195 -  0.00508)
#
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26182 ( 0.22104 -  0.31671)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73749 ( 0.68242 -  0.77833)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00070 ( 0.00027 -  0.00123)
#
# exp(N(mu,v)) has mean exp(mu+v/2) and variance exp(2(mu+v/2))(exp(v)-1)
mu1=jnp.concatenate([jnp.zeros(1),mu])
D1=jnp.concatenate([jnp.zeros([ncat,1]),jnp.concatenate([jnp.zeros([1,ncat-1]),D])],axis=1)
mu2=jnp.exp(mu1+jnp.diag(D1)/2)
D2=mu2[:,None]*mu2[None,:]*(jnp.exp(D1)-1)
t0=cpu()
print("Sampling exp(xx@newdata)")
oo=mvn.rvs(mean=mu2,cov=D2,size=numsamp)
print("Sampled")
pp=oo/oo.sum(axis=1)[:,None]
for i in range(ncat):
  print("Using %d samples of MVN: Prob(newdata = cat%d) = %8.5f (%8.5f - %8.5f)"%(numsamp,i,p[i],np.quantile(pp[:,i],(1-conf)/2),np.quantile(pp[:,i],(1+conf)/2)))
print("Time =",cpu()-t0)
print()

# Analytic CI pretending exp(xx@newdata) is normal and sum is positive
# df=2, mult=10
# Pretending denominator>0: Prob(newdata = cat0) =  0.26182 ( 0.15547 -  0.60668)
# Pretending denominator>0: Prob(newdata = cat1) =  0.73749 ( 0.39087 -  0.84382)
# Pretending denominator>0: Prob(newdata = cat2) =  0.00070 (-0.00194 -  0.00510)
#
# df=2, mult=100
# Pretending denominator>0: Prob(newdata = cat0) =  0.26182 ( 0.22104 -  0.31671)
# Pretending denominator>0: Prob(newdata = cat1) =  0.73749 ( 0.68243 -  0.77833)
# Pretending denominator>0: Prob(newdata = cat2) =  0.00070 ( 0.00027 -  0.00123)
#
zconf=norm.ppf((1+conf)/2)
mu2s=mu2.sum()
D2s=D2.sum(axis=1)
D2s2=D2s.sum()
a=mu2s**2-zconf**2*D2s2
for i in range(ncat):
  b=-2*mu2[i]*mu2s+2*zconf**2*D2s[i]
  c=mu2[i]**2-zconf**2*D2[i,i]
  di=b**2-4*a*c
  if a>0 and di>=0:
    s=sqrt(di)
    low=(-b-s)/(2*a)
    high=(-b+s)/(2*a)
    print("Pretending denominator>0: Prob(newdata = cat%d) = %8.5f (%8.5f - %8.5f)"%(i,p[i],low,high))
  else:
    print("Pretending denominator>0: Prob(newdata = cat%d) = %8.5f (n/a - n/a)"%(i,p[i]))
    

# MCMC over xx, i.e., exact multinomial posterior (with uniform prior)
# df: 2
# Multiplicity: 1
# Iterations: 24777382
# Acceptance ratio: 0.270102
# Time/iteration: 289.117us
# Category 0:   Median 0.193764      CI 0.007115 - 0.760380
# Category 1:   Median 0.804226      CI 0.235111 - 0.992809
# Category 2:   Median 0.000029      CI 0.000000 - 0.019679
#
# df: 2
# Multiplicity: 10
# Iterations: 70692757
# Acceptance ratio: 0.203377
# Time/iteration: 277.481us
# Category 0:   Median 0.255805      CI 0.132258 - 0.418642
# Category 1:   Median 0.743344      CI 0.580095 - 0.867277
# Category 2:   Median 0.000532      CI 0.000057 - 0.003591
#
# df: 2
# Multiplicity: 100
# Iterations: 24777382
# Acceptance ratio: 0.195965
# Time/iteration: 286.916us
# Category 0:   Median 0.261200      CI 0.216945 - 0.309441
# Category 1:   Median 0.738090      CI 0.689729 - 0.782447
# Category 2:   Median 0.000678      CI 0.000353 - 0.001262
#
it=ac=0
al=1
L=-NLL(xx)
t0=cpu()
l=[[] for i in range(ncat)]
nextprint=10000
while 1:
  xx1=mvn.rvs(mean=xx,cov=al*C)
  L1=-NLL1(xx1)
  if random()<exp(L1-L):
    xx=xx1;L=L1
    ac+=1
  yy=np.zeros([ncat,nvar])
  yy[1:,:]=xx.reshape([ncat-1,nvar])
  o=np.exp(yy@newdata)
  p=o/o.sum()
  for i in range(ncat): l[i].append(p[i])
  it+=1
  if it==nextprint:
    print("df:",df)
    print("Multiplicity:",mult)
    print("Iterations:",it)
    print("Acceptance ratio: %.6f"%(ac/it))
    print("Time/iteration: %.3fus"%((cpu()-t0)/it*1e6))
    n0=int((1-conf)/2*it)
    n1=int((1+conf)/2*it)
    med=it//2
    for i in range(ncat):
      l[i].sort()
      pp=np.array([l[i][med],l[i][n0],l[i][n1]])
      print("Category %d:"%i,"  Median %8.6f      CI %8.6f - %8.6f"%tuple(pp))
    print()
    nextprint=int(1.1*nextprint)
