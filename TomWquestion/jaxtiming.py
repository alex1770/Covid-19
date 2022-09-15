# https://jax.readthedocs.io/en/latest/notebooks/quickstart.html

from jax import jacfwd, jacrev
import jax.numpy as jnp
from jax import grad, jit, vmap
from jax import random
from math import log,exp
from scipy.stats import norm, multivariate_normal as mvn
from scipy.optimize import minimize
from time import process_time as cpu
#from jax.config import config
#config.update("jax_enable_x64", True)

# awk '{mpg=$1;cyl=$2;am=$9;lp=-8.3447988+0.3700268*mpg+0.7325165*(cyl==6)+0.7022486*(cyl==8);p=(am*exp(lp)+1-am)/(1+exp(lp));ll+=log(p);print am,p,lp}END{print ll}' < cardata

data=[]
with open('cardata') as fp:
  for x in fp:
    y=x.strip()
    if y=='' or y[0]=='#': continue
    z=y.split()
    # mpg (float), cyl (4, 6 or 8), am (0 or 1)
    data.append([float(z[0]),int(z[1]),int(z[8])])

data2=[];data3=[]
for (mpg,cyl,am) in data:
  data2.append([1,mpg,cyl==6,cyl==8])
  data3.append(am)# Using float doesn't make a huge difference
data2=jnp.array(data2)
data3=jnp.array(data3)
    
def NLL0(xx):
  LL=0
  for (mpg,cyl,am) in data:
    lp=xx[0]+xx[1]*mpg+xx[2]*(cyl==6)+xx[3]*(cyl==8)
    p=(am*jnp.exp(lp)+1-am)/(1+jnp.exp(lp))
    LL+=jnp.log(p)
  return -LL

def NLL(xx):
  p1=jnp.exp(data2@xx)
  return -(jnp.log((1-data3*(1-p1))/(1+p1))).sum()

d1=jit(grad(NLL))
d2=jit(jacfwd(jacrev(NLL)))

# [-8.3447988, 0.3700268, 0.7325165, 0.7022486]
#  -8.342127186627655   0.36992050557084866   0.7315879364536231   0.7016122038093382

bounds=[(-10,10), (0,1), (0,1), (0,1)]
t0=cpu()
res=minimize(NLL,[0.,0.5,0.5,0.5],jac=d1,bounds=bounds, method="SLSQP")
print(res.message)
print(cpu()-t0,"seconds")
print(res.fun)
print(res.x)
print()

t0=cpu()
xx=jnp.array([-8, 0.4, 0.7, 0.7])
for it in range(20):
  ee=jnp.linalg.solve(d2(xx),d1(xx))
  xx-=ee
  print(cpu()-t0,ee,xx)
  


