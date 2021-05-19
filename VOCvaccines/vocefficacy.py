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
fixr0=True#False

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
#VV=["B.1.1.7", "B.1.526"]#P.1"]#
#VV=["B.1.526", "B.1.1.7"]
VV=["B.1.1.7", "B.1.617.2"]

# Set up vax odds
oo=np.zeros([ng])
for g in range(ng):
  oo[g]=pv[g]/(1-pv[g])
  
# Set up incidence count matrix
numvar={}
nn=np.zeros([2,ng],dtype=int)
for (date,location,age,var) in zip(gg['Collection date'],gg['Location'],gg['Patient age'],gg['Pango lineage']):
  loc=location.split('/')
  if loc[1].strip()!=targcountry: continue
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
print(xx0,res.fun,res.nfev)

sd0=0.2
if fixr0:
  sd=[0,sd0,sd0]
else:
  sd=[sd0,sd0,sd0]

L0=LL(*xx0)

n=50
if fixr0:
  xx=np.array(xx0)
  for i1 in range(n):
    xx[1]=i1/n
    bounds=[(r0,r0),(xx[1],xx[1]),(0.001,1000)]
    res=minimize(NLL,xx,method="SLSQP",bounds=bounds,options={"maxiter":1000})
    if not res.success: raise RuntimeError(res.message)
    print("%7.4f   %7.2f"%(xx[1],res.fun+L0))
else:
  xx=np.array(xx0)
  for i0 in range(n):
    xx[0]=(i0+0.5)/n
    for i1 in range(n):
      xx[1]=(i1+0.5)/n
      bounds=[(xx[0],xx[0]),(xx[1],xx[1]),(0.001,1000)]
      res=minimize(NLL,xx,method="SLSQP",bounds=bounds,options={"maxiter":1000})
      if not res.success: raise RuntimeError(res.message)
      L=res.fun+L0
      if L<10:
        print(" %4.1f"%L,end='')
      else:
        print("   * ",end='')
    print()
