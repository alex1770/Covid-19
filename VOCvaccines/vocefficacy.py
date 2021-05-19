from stuff import *
import numpy as np
import sys
from math import log,exp
from scipy.optimize import minimize
from scipy.stats import norm
from random import random,seed
seed(42)

# Number of age bands
ng=7

targcountry="USA"
mindate="2021-03-15"
fixr0=True

# Vaccination data from https://covid.cdc.gov/covid-data-tracker/#vaccination-demographic
# As of 2021-05-16
# Todo: reduce to cover dates back to mindate
pv=[0.0399, 0.3659, 0.4492, 0.5148, 0.6190, 0.7966, 0.7759]

np.set_printoptions(precision=3,linewidth=120)

def agetogroup(x):
  try:
    a=int(x)
  except ValueError:
    return -1
  l=[18,30,40,50,65,75,999]
  for (i,y) in enumerate(l):
    if a<y: return i

gg=loadcsv("GISAID.ages.tsv",sep='\t')

# Variants to be compared
VV=["B.1.1.7", "B.1.526"]#P.1"]#
#VV=["B.1.526", "B.1.1.7"]
#VV=["B.1.1.7", "B.1.617.2"]

# Set up vax odds
oo=np.zeros([ng])
for g in range(ng):
  oo[g]=pv[g]/(1-pv[g])
  
# Set up incidence count matrix
numvar={}
nn=np.zeros([2,ng],dtype=int)
for (date,country,age,var) in zip(gg['date'],gg['country'],gg['age'],gg['pango_lineage']):
  if country!=targcountry: continue
  if date<mindate: continue
  g=agetogroup(age)
  if g<0: continue
  numvar[var]=numvar.get(var,0)+1
  if var not in VV: continue
  i=VV.index(var)
  nn[i,g]+=1

if 0:
  l=list(numvar)
  l.sort(key=lambda x:-numvar[x])
  for x in l:
    if numvar[x]<100: break
    print("%-11s %8d"%(x,numvar[x]))
  sys.exit(0)

print("Comparing variants",VV[0],"and",VV[1])
print()
print("Vaccination probabilities")
print(pv)
print()
print("Case counts")
print(nn)
print()

C=sum(nn[0])/sum(nn[1])# Abitrary constant that helps conditioning
def LL(r0,r1,D):
  DD=[C,D]
  rr=[r0,r1]
  tot=0
  for g in range(ng):
    l=[DD[i]*(1+oo[g]*rr[i]) for i in range(2)]
    s=sum(l)
    for i in range(2):
      tot+=nn[i,g]*log(l[i]/s)
  return tot

def NLL(xx): return -LL(*xx)

r0=0.1
xx=[r0,0.5,1]
if fixr0:
  bounds=[(r0,r0),(0.01,10),(0.001,1000)]
else:
  bounds=[(0.01,10),(0.01,10),(0.001,1000)]
res=minimize(NLL,xx,method="SLSQP",bounds=bounds,options={"maxiter":1000})
if not res.success: raise RuntimeError(res.message)
xx0=res.x
print(xx0,res.fun)

lcd=12
ncd=1<<lcd
hist=np.zeros(ncd)
cor=np.zeros([lcd+1,4])

sd0=0.2
if fixr0:
  sd=[0,sd0,sd0]
else:
  sd=[sd0,sd0,sd0]
nit=nacc=0
L0=LL(*xx0)
L=L0;xx=xx0

ln=0;n=0
while 1:
  dd=norm.rvs(size=3)
  yy=[x*exp(d*s) for (x,d,s) in zip(xx,dd,sd)]
  L2=LL(*yy)
  if random()<exp(L2-L): xx=yy;L=L2;nacc+=1
  n+=1
  if n==(1<<ln):
    r1=xx[1]
    cor[ln,0]+=1
    cor[ln,1]+=r1
    cor[ln,2]+=r1*r1
    ln+=1
    if ln>lcd:
      ln=0;n=0
      xx=xx0;L=L0
  nit+=1
  if n==0:
    print("%10d  %10d  %6.4f   %12g %12g   %12g"%(nit,nacc,nacc/nit,xx[0],xx[1],L-L0))
    for i in range(lcd+1):
      cc=cor[i]/cor[i,0]
      v1=cc[2]-cc[1]**2
      v2=(cor[i,2]-cor[i,1]**2/cor[i,0])/(cor[i,0]-1)
      print("%3d  %12g  %12g"%(i,v1,v2))
    print()
while 1:
  dd=norm.rvs(size=3)
  yy=[x*exp(d*s) for (x,d,s) in zip(xx,dd,sd)]
  L2=LL(*yy)
  if random()<exp(L2-L): xx=yy;L=L2;nacc+=1
  if nit>=ncd:
    r1=xx[1]
    for i in range(lcd+1):
      q1=hist[(nit-(1<<i))&(ncd-1)]
      cor[i,0]+=1
      cor[i,1]+=r1
      cor[i,2]+=r1*r1
      cor[i,3]+=q1*r1
    hist[nit&(ncd-1)]=r1
  nit+=1
  if nit%1000==0:
    print("%10d  %6.4f   %12g %12g   %12g"%(nacc,nacc/nit,xx[0],xx[1],L-L0))
    for i in range(lcd+1):
      cc=cor[i]/cor[i,0]
      v1=cc[2]-cc[1]**2
      v2=cc[3]-cc[1]**2
      print("%3d  %12g  %12g  %12g"%(i,v1,v2,v2/v1))


"""
fixr0=True
  13844776  0.1331            0.1    0.0501897       -2.20353
  0   0.000280851   0.000276936      0.986062
  1   0.000280851   0.000273084      0.972346
  2   0.000280851   0.000265566      0.945576
  3   0.000280851   0.000251225      0.894513
  4   0.000280851    0.00022511       0.80153
  5   0.000280851   0.000181637      0.646737
  6   0.000280851   0.000120131      0.427741
  7   0.000280851   5.50729e-05      0.196093
  8   0.000280851   1.29277e-05     0.0460306
  9   0.000280851   5.29983e-07    0.00188706
 10   0.000280851  -6.55516e-08  -0.000233403
 11   0.000280851   1.89097e-07     0.0006733
 12   0.000280851  -1.04112e-07  -0.000370702
KeyboardInterrupt
>>> nit
103990271
>>> 
"""
