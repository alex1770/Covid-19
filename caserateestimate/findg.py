import requests,sys
import numpy as np
from math import log,exp
from stuff import *

mindate="2022-04-01"
if len(sys.argv)>1: mindate=sys.argv[1]
print("mindate =",mindate,file=sys.stderr)

population=55.98e6

now=Date(datetime.datetime.now().strftime("%Y-%m-%d"))

data0=getcases_raw("2022-05-25",location="England")
data=getcases_raw("2022-06-01",location="England")

data0={Date(d):c for (d,c) in data0.items() if d>=mindate}
data={Date(d):c for (d,c) in data.items() if d>=mindate}
last=max(data)

# Correct recent incomplete entries, based on previous week
for d in Daterange(last-6,last+1):
  data[d]=round(data[d]*data[d-7]/data0[d-7])

date0=min(data)
# Check dates are continguous
assert [x-date0 for x in list(data)]==list(range(len(data)))

counts0=list(data.values())
counts=weekdayadj(counts0)

# Discard last entry as unreliable
counts=counts[:-1]
n=len(counts)
extrap=21

exps=[
  ("BA.2",      0,     Date("2022-01-01")),
  ("BA.2.12.1", 0.111, int(Date("2022-06-08"))+0.1),
  ("BA.4",      0.114, int(Date("2022-06-07"))+0.6),
  ("BA.5",      0.133, int(Date("2022-06-04"))+0.3)
]

# Get variant factors vs BA.2 baseline
#  BA.2   ==BA.2.12.1==    ====BA.4=====    ====BA.5=====
#     1 + exp(r1*(t-t1)) + exp(r2*(t-t2)) + exp(r3*(t-t3))
vf=[]
for i in range(n+extrap):
  d=date0+i
  s=0
  for (name,g,cross) in exps: s+=exp(g*(int(d)-cross))
  vf.append(s)

for (source,fn) in [(counts0,"residnonadj"), (counts,"resid")]:
  with open(fn,'w') as fp:
    for i in range(n):
      d=date0+i
      c=source[i]
      s=vf[i]
      print(d,"%7d %12g %7.3f"%(c,s,log(c/s)),file=fp)
# gnuplot makeresidgraph.gpl
# gnuplot makeresidgraphnonadj.gpl

g0values=[-0.065,-0.053]
intercept=(Date("2022-05-29"),7.75)
with open("extrap","w") as fp:
  print("#     Date       Actual       Tot_l       Tot_h",end="",file=fp)
  for (name,g,cross) in exps:
    print("%12s%12s"%(name+"_l",name+"_h"),end="",file=fp)
  print(file=fp)
  for i in range(n+extrap):
    d=date0+i
    if i<n: c="%12g"%counts[i]
    else: c="%12s"%"-"
    print(d,c,end="",file=fp)
    for g0 in g0values:
      est=vf[i]*exp(intercept[1]+g0*(d-intercept[0]))
      print("%12g"%est,end="",file=fp)
    for (name,g,cross) in exps:
      for g0 in g0values:
        est=exp(g*(int(d)-cross))*exp(intercept[1]+g0*(d-intercept[0]))
        print("%12g"%est,end="",file=fp)
    print(file=fp)
# gnuplot makeextrapgraphs.gpl

print("     Variant    Daily growth     Doub/Half-time       Approx R (*)   (* using gen time = 4 days)")
print("     =======  ==============     ==============      =============")
for (name,g,cross) in exps:
  print("%12s"%name,end="")
  for g0 in g0values:
    print("  %6.3f"%(g0+g),end="")
  print("   ",end="")
  for g0 in g0values:
    print("  %6.1f"%(log(2)/(g0+g)),end="")
  print("   ",end="")
  for g0 in g0values:
    print("  %6.2f"%(exp(4*(g0+g))),end="")
  print()
  
