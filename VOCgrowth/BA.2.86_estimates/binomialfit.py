# Binomial fit, with likelihood output of growth rate, and possibly introduction day
# What happens if you add a day-of-introduction parameter independent of the constant offset?
# Explain (a_i, b_i)   a_i=non-variant, b_i=variant,
# by assuming a_i+b_i is given and that b_i ~ B(a_i+b_i, p_i)
# where p_i/(1-p_i) = exp(c + (i-i0)*g)  if i>=i0
#                   = 0                  if i<i0
# But c isn't allowed to be arbitrary, because on the day of introduction
# we expect the chance of seeing the variant to be on the order of (a few)/number_of_cases, i.e.,
# c should be about -log(number_of_cases/(a few)).
# So in that case, we don't need to add the 0 for i<i0, because p_i will be
# essentially zero for i<i0 anyway (assuming g>=0), given that we're sequencing
# far fewer than the number of infections. So we're back to normal binomial regression.
#
# So here we do normal binomial regression, where p_i/(1-p_i) = exp(c + i*g)
#
# Rather than integrating over the nuisance parameter, c, uniformly, we'd like to get a bit more power by using a suitable prior.
# Assume 1 virus on the day of introduction, so c+introday*g = -log(ipd), where ipd=infections per day.
# (This assumes that we're comparing variant with everything else. Need to adjust if comparing variant to variant.)
# Don't know ipd, but a reasonable range for the territories we are using (populous countries India, USA, China are broken down into states)
# is 3000 - 600000, so treat log(ipd) as U[log(3000),log(600000)].
# The prior for the introduction of BA.2.86 is [2023-03-01, variantfirstseen]. (Would need to modify for other variants.)
# So, for a given g, we assume prior for c is -(U[2023-03-01, variantfirstseen] + U[log(3000),log(600000)]).

import numpy as np
from math import exp,log,sqrt,floor
from stuff import *
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-m', '--maxg',     type=float,default=0.2,     help="Maximum daily logarithmic growth rate considered (effectively the prior is U[0,maxg])")
parser.add_argument('-f', '--minintrodate',  default="2023-03-01",  help="Earliest possible introduction date of variant")
parser.add_argument('-w', '--writegraph', action="store_true",      help="Whether to write graph output file")
parser.add_argument('countfilenames',   nargs='*',                  help="Name of file containing counts of non-variant, variant")
args=parser.parse_args()

# Odds of variant:non-variant are modelled as exp(c+i*g), i=day number
# p_i/(1-p_i) = exp(c + i*g), i=0, 1, ..., ndays-1
def LL(g,c,V0,V1):
  ndays=len(V0)
  logodds=c+np.arange(ndays)*g
  Z=np.log(1+np.exp(logodds))
  return np.sum(V1*(logodds-Z)-V0*Z)

dg=0.001
def g2bin(g): return int(floor(g/dg+0.5))
def bin2g(b): return b*dg

minintrodate=datetoday(args.minintrodate)
minipd=3000
maxipd=600000
dc=0.25

def getlik(countfile):

  source="GISAID"
  # Variant0, Variant1 counts by day
  N0=[];N1=[];DT=[]
  with open(countfile) as fp:
    for x in fp:
      if x[0]=='#': source=x[1:].strip();continue
      y=x.strip().split()
      d=datetoday(y[0])
      v0,v1=float(y[1]),float(y[2])
      N0.append(v0);N1.append(v1);DT.append(d)
  
  minday=min(DT)
  maxday=max(DT)+1
  ndays=maxday-minday
  V0=np.zeros(ndays)
  V1=np.zeros(ndays)
  for (v0,v1,dt) in zip(N0,N1,DT):
    V0[dt-minday]=v0
    V1[dt-minday]=v1
  
  firstseen=min(i for i in range(ndays) if V1[i]>0)

  l=[]
  ming,maxg=0,args.maxg
  mingq,maxgq=bin2g(g2bin(ming)),bin2g(g2bin(maxg))
  logsum={};logmax={}
  # Take expectation over c ~ -U[minintrodate, firstseen]*g - U[log(minipd),log(maxipd)]
  A1=log(minipd);B1=log(maxipd)
  for g in np.arange(mingq,maxgq,dg):
    A0=g*(minintrodate-minday);B0=g*firstseen
    if A0>B0: raise RuntimeError("First seen variant before min intro date in "+countfile)
    # c ~ -(U[A0,B0]+U[A1,B1])
    ll=[]
    totw=0
    for mc in np.arange(A0+A1,B0+B1,dc):
      weight=min(mc-(A0+A1),B0+B1-mc,B0-A0,B1-A1)+1e-9
      totw+=weight
      ll.append(LL(g,-mc,V0,V1)+log(weight))
    logmax[g]=max(ll)-log(totw)
    logsum[g]=log(sum(exp(x-logmax[g]) for x in ll))+logmax[g]-log(totw)

  return source,logsum,logmax

def printstats(source,ll,desc):
  print(f"Relative growth rate estimate of BA.2.86 at {UKdatetime()[0]}. Data from {source}.")
  mx=max(ll.values())
  el={g:exp(ll[g]-mx) for g in ll}
  tot=sum(el.values())
  if args.writegraph:
    with open("outlik_"+desc,"w") as fp:
      for g in sorted(list(el)):
        print("%10.6f %12g"%(g,el[g]/tot/dg),file=fp)
  thr=[0.05,0.1,0.25,0.5,0.9,0.95,0.99]
  i=0;t=0;s=0
  for g in sorted(list(el)):
    p=el[g]/tot
    t+=p
    s+=p*g
    while i<len(thr) and t>thr[i]:
      print("%s: %2.0f%% point at daily logarithmic growth %5.3f = %+8.0f%% per week"%(desc,thr[i]*100,g,(exp(g*7)-1)*100))
      i+=1
  print("%s: mean      at daily logarithmic growth %5.3f = %+8.0f%% per week"%(desc,s,(exp(s*7)-1)*100))
  print()
    

post={}
nn={}
sources=set()
for countfile in args.countfilenames:
  source,logsum,logmax=getlik(countfile)
  printstats(source,logsum,countfile)
  sources.add(source)
  #print("%12g %12g %12g %12g %12g"%(logsum[0],logsum[0.05],logsum[0.1],logsum[0.15],logsum[0.199]))
  for g in logsum:
    post[g]=post.get(g,0)+logsum[g]
    nn[g]=nn.get(g,0)+1

n=len(args.countfilenames)
post2={g:post[g] for g in post if nn[g]==n}
printstats(" and ".join(sources),post2,"Combined")
