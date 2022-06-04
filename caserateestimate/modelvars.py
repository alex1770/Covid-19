import requests,sys
import numpy as np
from math import log,exp
from stuff import *

# Project back estimated infection count onto a series of variants with hypothesised growth relationships between them

population=55.98e6

inf={}
with open('estdailyinfectionsEngland') as fp:
  for x in fp:
    y=x.split()
    inf[Date(y[0])]=float(y[1])

date0=min(inf)
counts=list(inf.values())
n=len(counts)

vars=["BA.1","BA.1.1","BA.2","BA.2.12.1","BA.4","BA.5"]

# Relative growths
# Entry name0, name1, g, t0 in relexps means we expect name1/name0 = exp(g*(t-t0))
relexps=[
  ("BA.1", "BA.1.1",    0.039, int(Date("2022-02-05"))+0.4),
  ("BA.1", "BA.2",      0.119, int(Date("2022-02-12"))+0.0),
  ("BA.2", "BA.2.12.1", 0.111, int(Date("2022-06-08"))+0.1),
  ("BA.2", "BA.4",      0.114, int(Date("2022-06-07"))+0.6),
  ("BA.2", "BA.5",      0.133, int(Date("2022-06-04"))+0.3)
]

# Entry name:(g,h) in exps means we expect name/baseline = exp(g*t+h)
baseline="BA.1"
exps={baseline: (0,0)}
for (name0,name1,g,t0) in relexps:
  exps[name1]=(exps[name0][0]+g,exps[name0][1]-g*t0)
  

# Get variant factors vs baseline
vf=[]
for i in range(n):
  d=date0+i
  s=0
  for (g,h) in exps.values(): s+=exp(g*int(d)+h)
  vf.append(s)

with open("variantest",'w') as fp:
  print("#Date      Infections",end="",file=fp)
  for name in exps:
    print("  %10s"%name,end="",file=fp)
  print(file=fp)
  for i in range(n):
    d=date0+i
    c=counts[i]
    s=vf[i]
    print(d,"%10d"%c,end="",file=fp)
    for (g,h) in exps.values():
      x=exp(g*int(d)+h)
      print("  %10.1f"%(x/s*c),end="",file=fp)
    print(file=fp)

print("     Variant  Relative growth")
print("     =======  ===============")
for name in exps:
  print("%12s           %6.3f"%(name,exps[name][0]))
