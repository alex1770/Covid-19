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
# essentially zero for i<i0 anyway (assuming g>=0 before it was first seen), given that we're sequencing
# far fewer than the number of infections. So we're back to normal binomial regression.
#
# So here we do normal binomial regression, where p_i/(1-p_i) = exp(c + i*g)
#
# Rather than integrating over the nuisance parameter, c, uniformly, we'd like to get a bit more power by using a suitable prior.
# Assume 1 virus on the day of introduction, so at p_{introday}=1/(ipd+1), where ipd=infections per day of non-variant, so
# at i=introday, 1/ipd = p_i/(1-p_i) = exp(c + i*g), so c+introday*g = -log(ipd), 
# (This assumes that we're comparing a variant with everything else. If comparing variant to variant, then ipd could get much lower, so it's
# best to organise it so that the second variant is the rarer one, and minipd is sufficiently lower to account for min first variant numbers.)
# Don't know ipd, but a reasonable range for the territories we are using (populous countries India, USA, China are broken down into states)
# is 3000 - 600000, so could treat log(ipd) as U[log(3000),log(600000)].
# The prior for the introduction of BA.2.86 is [2023-03-01, variantfirstseen]. (Would need to modify for other variants.)
# So, for a given g, we'd assume prior for c is -(U[2023-03-01, variantfirstseen] + U[log(3000),log(600000)]).
# However, it's not sensitive to the lower end of ipd (3000 in the above), so we can make it small enough to accommodate other situations,
# like variant vs variant.


import numpy as np
from math import exp,log,sqrt,floor
from stuff import *
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-l', '--ming',     type=float,default=-0.05,   help="Minimum daily logarithmic growth rate considered")
parser.add_argument('-m', '--maxg',     type=float,default=0.15,    help="Maximum daily logarithmic growth rate considered")
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

dg=0.0005
def g2bin(g): return int(floor(g/dg+0.5))
def bin2g(b): return b*dg

minintrodate=datetoday(args.minintrodate)
minipd=1      # min and max infections per day of non-variant - see above
maxipd=600000 #
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

  if 0:#alter
    for i in range(len(N1)):
      if N1[i]>0: break
    else: assert 0
    i+=1
    DT=DT[i:];N0=N0[i:];N1=N1[i:]
  if sum(N1)==0: return (source,)+({g:1 for g in np.arange(args.ming,args.maxg,dg)},)*2
  
  minday=min(DT)
  maxday=max(DT)+1
  ndays=maxday-minday
  V0=np.zeros(ndays)
  V1=np.zeros(ndays)
  for (v0,v1,dt) in zip(N0,N1,DT):
    V0[dt-minday]=v0
    V1[dt-minday]=v1
  #if V0.sum()<V1.sum(): print("Warning: baseline variant should not be smaller than new variant in input",countfile,file=sys.stderr)# alter
  
  firstseen=min(i for i in range(ndays) if V1[i]>0)

  l=[]
  ming,maxg=args.ming,args.maxg
  mingq,maxgq=bin2g(g2bin(ming)),bin2g(g2bin(maxg))
  logsum={};logmax={}
  # Take expectation over c ~ -U[minintrodate, firstseen]*g - U[log(minipd),log(maxipd)]
  A1=log(minipd);B1=log(maxipd)
  if minintrodate-minday>firstseen: raise RuntimeError("First seen variant before min intro date in "+countfile)
  n=int(floor((B1-A1+(firstseen-(minintrodate-minday))*max(abs(mingq),abs(maxgq)))/dc))+1
  for g in np.arange(mingq,maxgq,dg):
    A0=g*(minintrodate-minday);B0=g*firstseen
    if g<0: A0,B0=B0,A0
    # c ~ -(U[A0,B0]+U[A1,B1])
    ll=[]
    totw=0
    for i in range(n):
      mc=(B0+B1-A0-A1)*(i+0.5)/n+(A0+A1)
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
    dn=os.path.dirname("outlik_"+desc)
    if dn!="": os.makedirs(dn,exist_ok=True)
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
      print("%s: %2.0f%% point at daily logarithmic growth %6.3f = %+8.0f%% per week"%(desc,thr[i]*100,g,(exp(g*7)-1)*100))
      i+=1
  print("%s: mean      at daily logarithmic growth %6.3f = %+8.0f%% per week"%(desc,s,(exp(s*7)-1)*100))
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
