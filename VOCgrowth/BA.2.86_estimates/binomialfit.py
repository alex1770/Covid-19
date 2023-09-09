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
# So here we do normal binomial regression,
# where p_i/(1-p_i) = exp(c + i*g)
#
# outputting a likelihood for g, assuming g>=0.

import numpy as np
from math import exp,log,sqrt,floor
from stuff import *
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-i', '--ipd',      type=float,default=10000,   help="Guessed total number of infections per day (shouldn't matter much)")
parser.add_argument('-m', '--maxg',     type=float,default=0.2,     help="Maximum daily logarthmic growth rate considered (effectively the prior is U[0,maxg])")
parser.add_argument('-f', '--minintrodate',  default="2019-01-01",  help="Earliest possible introduction date of variant")
parser.add_argument('-w', '--writegraph', action="store_true",      help="Whether to write graph output file")
parser.add_argument('countfilenames',   nargs='*',                  help="Name of file containing counts of non-variant, variant")
args=parser.parse_args()

# Odds of variant:non-variant are modelled as exp(c+i*g), i=day number
# Assuming on day of introduction odds are 1/ipd
# So c + intro*g = -log(ipd)
# p_i/(1-p_i) = exp(c + i*g), i=0, 1, ..., ndays-1
def LL(g,intro,ipd,V0,V1):
  c=-log(ipd)-intro*g
  ndays=len(V0)
  logodds=c+np.arange(ndays)*g
  Z=np.log(1+np.exp(logodds))
  return np.sum(V1*(logodds-Z)-V0*Z)

dg=0.001
def g2bin(g): return int(floor(g/dg+0.5))
def bin2g(b): return b*dg

minintrodate=datetoday(args.minintrodate)

def getlik(countfile):

  # Variant0, Variant1 counts by day
  N0=[];N1=[];DT=[]
  with open(countfile) as fp:
    for x in fp:
      if x[0]=='#': continue
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
  dsum={}# Integrating over nuisance parameter (c or intro)
  dmax={}# Maximising over nuisance parameter (c or intro)
  for intro in np.arange(max(0,minintrodate-minday),firstseen,0.25):
    ming,maxg = 0,args.maxg
    for g in np.arange(bin2g(g2bin(ming)),bin2g(g2bin(maxg)),dg):
      ll=LL(g,intro,args.ipd,V0,V1);el=exp(ll)
      dsum[g]=dsum.get(g,0)+el
      dmax[g]=max(dmax.get(g,0),el)
      #print("%8.3f %10.6f %10.6f %12g"%(intro,g,ll,el),file=fp)

  return dsum,dmax

sources=set()
def printstats(d,desc):
  if desc=="Combined":
    source=" and ".join(sources)
  else:
    source="COG-UK" if "COG" in desc else "GISAID"
    sources.add(source)
  print(f"Relative growth rate estimate of BA.2.86 at {UKdatetime()[0]}. Data from {source}.")
  tot=sum(d.values())
  if args.writegraph:
    with open("outlik_"+desc,"w") as fp:
      for g in sorted(list(d)):
        print("%10.6f %12g"%(g,d[g]/tot/dg),file=fp)
  thr=[0.05,0.1,0.25,0.5,0.9,0.95,0.99]
  i=0;t=0;s=0
  for g in sorted(list(d)):
    p=d[g]/tot
    t+=p
    s+=p*g
    while i<len(thr) and t>thr[i]:
      print("%s: %2.0f%% point at daily logarithmic growth %5.3f = %+8.0f%% per week"%(desc,thr[i]*100,g,(exp(g*7)-1)*100))
      i+=1
  print("%s: mean      at daily logarithmic growth %5.3f = %+8.0f%% per week"%(desc,s,(exp(s*7)-1)*100))
  print()
    

post={}
nn={}
for countfile in args.countfilenames:
  dsum,dmax=getlik(countfile)
  printstats(dsum,countfile)
  for g in dsum:
    post[g]=post.get(g,1)*dsum[g]
    nn[g]=nn.get(g,0)+1

n=len(args.countfilenames)
post2={g:post[g] for g in post if nn[g]==n}
printstats(post2,"Combined")
