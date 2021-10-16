import sys
from stuff import *
import numpy as np
from math import log,sqrt

mindate='2021-07-01'
maxdate='9999-12-31'

minday=datetoday(mindate)

#fp=sys.stdin
fp=open('cog_metadata.csv','r')#alter
mutcounts={}
daycounts={}
mutdaycounts={}
for (dt,lin,var,mut) in csvrows_it(fp,['sample_date','lineage','scorpio_call','mutations']):
  if dt>maxdate: continue
  day=datetoday(dt)
  #if day<minday-30: break# Can't use this cutoff in a simple way because of lots of erroneous dates
  if day<minday: continue
  daycounts[day]=daycounts.get(day,0)+1
  muts=mut.split('|')
  for mut in muts:
    mutcounts[mut]=mutcounts.get(mut,0)+1
    if mut not in mutdaycounts: mutdaycounts[mut]={}
    mutdaycounts[mut][day]=mutdaycounts[mut].get(day,0)+1
    
growth={}
for mut in mutcounts:
  if mutcounts[mut]<100: continue
  #if mut!='synSNP:C27143T': continue
  m=np.zeros([2,2])
  r=np.zeros(2)
  for day in mutdaycounts[mut]:
    x=day-minday
    vm=mutdaycounts[mut][day]
    v0=daycounts[day]
    if v0>0 and vm>0:
      y=log(vm/v0)
      w=1/(1/v0+1/vm)
      m[0,0]+=w
      m[0,1]+=w*x
      m[1,0]+=w*x
      m[1,1]+=w*x*x
      r[0]+=w*y
      r[1]+=w*x*y
  mi=np.linalg.pinv(m)
  c=np.linalg.solve(m,r)
  cv=[mi[0,0],mi[1,1]]# These should be the variances of c[0],c[1]
  #print(c[1],sqrt(cv[1]))
  #print("Continuous growth rate = %.4f/day"%(c[1]))
  #print("Crossover on",daytodate(datetoday(dt[0])+int(-c[0]/c[1]+.5)))
  #print("%-20s  %7.4f   %s"%(mut,c[1],daytodate(minday+int(-c[0]/c[1]+.5))))
  growth[mut]=c[1],sqrt(cv[1])

sd=6
l=[mut for mut in growth if growth[mut][1]>0]
#l.sort(key=lambda x:-(growth[x][0]-sd*growth[x][1]))
l.sort(key=lambda x:-(growth[x][0]/growth[x][1]))
nm=0
for mut in l:
  gr=growth[mut]
  (g,gl,gh)=(gr[0],gr[0]-sd*gr[1],gr[0]+sd*gr[1])
  if gl<0: break
  print("%-20s  %6.3f (%6.3f - %6.3f) %6.2f   %8d"%(mut,g,gl,gh,gr[0]/gr[1],mutcounts[mut]))
  nm+=1

nmd=min(nm,10)
print()
days=set()
for mut in l[:nmd]: days.update(mutdaycounts[mut])
print("                All",end='')
for mut in l[:nmd]: print(" %20s"%mut,end='')
print()
days=sorted(days)
for day in days:
  v0=daycounts[day]
  print(daytodate(day),"  %6d"%v0,end='')
  for mut in l[:nmd]:
    vm=mutdaycounts[mut].get(day,0)
    if vm>0 and v0>0: print("         %6.3f"%(log(vm/v0)),end='')
    else: print("              -",end='')
    print(" %5d"%vm,end='')
  print()

def grad(ll):
  m=np.zeros([2,2])
  r=np.zeros(2)
  for (x,v0,vm) in ll:
    if v0>0 and vm>0:
      y=log(vm/v0)
      w=1/(1/v0+1/vm)
      m[0,0]+=w
      m[0,1]+=w*x
      m[1,0]+=w*x
      m[1,1]+=w*x*x
      r[0]+=w*y
      r[1]+=w*x*y
  c=np.linalg.solve(m,r)
  mi=np.linalg.pinv(m)
  c=np.linalg.solve(m,r)
  return c[1],mi[1,1]

from scipy.stats import poisson,gamma
def perturb(ll):
  mm=[]
  for (x,v0,vm) in ll:
    #v0p=poisson.rvs(v0)
    #vmp=poisson.rvs(vm)
    v0p=gamma.rvs(v0) if v0>0 else 0
    vmp=gamma.rvs(vm) if vm>0 else 0
    mm.append((x,v0p,vmp))
  return mm

from random import randrange
for tr in range(100):
  mut=l[randrange(nm)] if tr>0 else 'synSNP:A4060T'
  ll=[]
  for day in days:
    v0=daycounts[day]
    vm=mutdaycounts[mut].get(day,0)
    #print(daytodate(day),"  %6d   %6d"%(v0,vm))
    ll.append((day-days[0],v0,vm))
  
  mu,va=grad(ll)
  sd0=sqrt(va)
  #print(mut,mu,sd0)
  nit=250
  s0=s1=s2=s3=0
  for it in range(nit):
    mm=perturb(ll)
    x=grad(mm)[0]
    s0+=1
    s1+=x
    s2+=x*x
    s3+=(x-mu)**2
  sd1=sqrt(s3/s0)
  sd2=sqrt((s2-s1**2/s0)/(s0-1))
  print("%-20s %9.7f | %9.7f %9.7f %9.7f | %9.6f %9.6f"%(mut,mu,sd0,sd1,sd2,sd1/sd0,sd2/sd0))

