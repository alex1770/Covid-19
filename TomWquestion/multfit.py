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

# Modification of cardata to keep MLE bounded.
# Otherwise we have problems because (e.g.,) mpg<21.4 ==> cyl=4, mpg>21.4 ==> cyl>4
# So you can add an arbitrary multiple of (mpg-21.4) to the predictors for cyl=6,8 and always improve the likelihood
# which means there is no MLE: the betas are unbounded.
suffix="_bounded"

data=[]
with open('cardata'+suffix) as fp:
  for x in fp:
    y=x.strip()
    if y=='' or y[0]=='#': continue
    z=y.split()
    # mpg (float), cyl (4, 6 or 8), am (0 or 1)
    data.append([float(z[0]),int(z[1]),int(z[8])])

nsmpg2=[]
with open('nsmpg2'+suffix) as fp:
  for x in fp:
    nsmpg2.append([float(z) for z in x.strip().split()])
    
nsmpg1=[]
with open('nsmpg1'+suffix) as fp:
  for x in fp:
    nsmpg1.append([float(z) for z in x.strip().split()])

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
# NLL = 11.644849119971521*mult
# [[  0.         0.         0.         0.      ]
#  [-26.976948  20.418986 -84.07624    0.028558]
#  [ 14.021035 -27.081224  -3.297667  -5.352148]]

newdata=data_in[0]
print(yy@newdata)
# [ 0.        1.033937 -5.909332]
o=np.exp(yy@newdata)
p=o/o.sum()
print("Newdata probs =",p)
numsamp=1000000

if 0:
  # See immediately below this for more efficient version of the same thing
  # CI with MVN (on xx, the log scale)
  # df=2, mult=1
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
# Using 500000000 samples of MVN: Prob(newdata = cat0) = 0.262135 (0.030118 - 0.781554)
# Using 500000000 samples of MVN: Prob(newdata = cat1) = 0.737154 (0.176761 - 0.967048)
# Using 500000000 samples of MVN: Prob(newdata = cat2) = 0.000711 (0.000001 - 0.254774)
#
# df=2, mult=10
# Using 500000000 samples of MVN: Prob(newdata = cat0) = 0.262135 (0.142799 - 0.430576)
# Using 500000000 samples of MVN: Prob(newdata = cat1) = 0.737154 (0.567646 - 0.856529)
# Using 500000000 samples of MVN: Prob(newdata = cat2) = 0.000711 (0.000094 - 0.005171)
#
# df=2, mult=100
# Using 500000000 samples of MVN: Prob(newdata = cat0) = 0.262135 (0.218520 - 0.310958)
# Using 500000000 samples of MVN: Prob(newdata = cat1) = 0.737154 (0.688176 - 0.780843)
# Using 500000000 samples of MVN: Prob(newdata = cat2) = 0.000711 (0.000377 - 0.001338)

mu=xx.reshape([ncat-1,nvar])@newdata
D=(C.reshape([ncat-1,nvar,ncat-1,nvar])*newdata[None,:,None,None]*newdata[None,None,None,:]).sum(axis=[1,3])
t0=cpu()
print("Sampling xx@newdata")
samp=mvn.rvs(mean=mu,cov=D,size=numsamp)
print("Sampled")
oo=jnp.concatenate([jnp.ones([numsamp,1]),jnp.exp(samp)],axis=1)
pp=oo/oo.sum(axis=1)[:,None]
for i in range(ncat):
  print("Using %d samples of MVN: Prob(newdata = cat%d) = %8.6f (%8.6f - %8.6f)"%(numsamp,i,p[i],np.quantile(pp[:,i],(1-conf)/2),np.quantile(pp[:,i],(1+conf)/2)))
print("Time =",cpu()-t0)
print()

# (Consider also doing this by direct integral: either in normal approx, or fully with multinomial)

# CI with MVN on exp(xx@newdata) (kind of stupid, but including anyway because can do it semi-analytically)
# df=2, mult=1
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26213 (-0.29085 -  0.29333)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73715 (-2.79183 -  2.90323)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00071 (-2.02956 -  3.91552)

# df=2, mult=10
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26213 ( 0.15554 -  0.60228)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73715 ( 0.39525 -  0.84373)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00071 (-0.00198 -  0.00517)
#
# df=2, mult=100
# Using 100000000 samples of MVN: Prob(newdata = cat0) =  0.26213 ( 0.22137 -  0.31697)
# Using 100000000 samples of MVN: Prob(newdata = cat1) =  0.73715 ( 0.68215 -  0.77799)
# Using 100000000 samples of MVN: Prob(newdata = cat2) =  0.00071 ( 0.00027 -  0.00125)
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
# Pretending denominator>0: Prob(newdata = cat0) =  0.26213 ( 0.15579 -  0.60593)
# Pretending denominator>0: Prob(newdata = cat1) =  0.73715 ( 0.39158 -  0.84348)
# Pretending denominator>0: Prob(newdata = cat2) =  0.00071 (-0.00197 -  0.00519)
#
# df=2, mult=100
# Pretending denominator>0: Prob(newdata = cat0) =  0.26213 ( 0.22137 -  0.31698)
# Pretending denominator>0: Prob(newdata = cat1) =  0.73715 ( 0.68213 -  0.77799)
# Pretending denominator>0: Prob(newdata = cat2) =  0.00071 ( 0.00027 -  0.00125)
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
#
# df: 2
# Multiplicity: 1
# Iterations: 36276564
# Acceptance ratio: 0.269666
# Time/iteration: 285.536us
# Category 0:   Median 0.194395      CI 0.007197 - 0.759243
# Category 1:   Median 0.803560      CI 0.236196 - 0.992717
# Category 2:   Median 0.000030      CI 0.000000 - 0.020185
# 
# df: 2
# Multiplicity: 10
# Iterations: 36276564
# Acceptance ratio: 0.203362
# Time/iteration: 286.034us
# Category 0:   Median 0.256059      CI 0.132571 - 0.418683
# Category 1:   Median 0.743074      CI 0.580025 - 0.866950
# Category 2:   Median 0.000545      CI 0.000058 - 0.003662
# 
# df: 2
# Multiplicity: 100
# Iterations: 36276564
# Acceptance ratio: 0.195855
# Time/iteration: 286.830us
# Category 0:   Median 0.261572      CI 0.217324 - 0.309741
# Category 1:   Median 0.737700      CI 0.689417 - 0.782064
# Category 2:   Median 0.000693      CI 0.000361 - 0.001288

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
