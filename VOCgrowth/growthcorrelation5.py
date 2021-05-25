from stuff import *
import sys
from scipy.optimize import minimize
from math import log,sqrt
from scipy.stats import norm

# Estimate growth advantage of B.1.617.2 over B.1.1.7 by correlating, over LTLAs, (a) the growth over
# a pair of weeks and (b) the relative prevalence of B.1.617.2 as estimated by Sanger sequencing data.

# Get ltla.csv from https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv
# Get Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv

apicases=loadcsv("ltla.csv")
sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')

lastsanger=max(sanger['WeekEndDate'])
#sangerday=datetoday(lastsanger)-7# Use penultimate week, because last week is incomplete
sangerday=datetoday(lastsanger)
sangerdate=daytodate(sangerday)
date0='2021-05-04'
date1=daytodate(datetoday(max(apicases['date']))-1)# Require dates to be less than date1, so this ignores last 2 days of cases-by-specimen-day

mgt=5# Mean generation time in days
exclude=set()#set(["E08000001"])

day0=datetoday(date0)
day1=datetoday(date1)
cases={}
for (ltla,date,n) in zip(apicases['areaCode'],apicases['date'],apicases['newCasesBySpecimenDate']):
  if ltla not in cases: cases[ltla]=[0]*(day1-day0+13)
  day=datetoday(date)
  if day>=day0-13 and day<day1: cases[ltla][day-day0+13]+=n

ltlas=set()
sangnum={}
for (date,ltla,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
  if ltla in cases and ltla not in exclude and date==sangerdate:
    if ltla not in sangnum: sangnum[ltla]=[0,0]
    if var=="B.1.617.2": sangnum[ltla][0]+=n
    sangnum[ltla][1]+=n
ltlas=sorted(list(sangnum))

# Log likelihood of reproduction numbers Q and R.
# Given r B.1.617.2 sequences found out of n,
# and given an overall increase in cases from A at time t1 to B at time t2 (times relative to Sanger sample),
# find likelihood of most likely p,a,b, where r~B(n,p), A~Po(a), B~Po(b), subject to ((1-p)Q^t2+pR^t2)/((1-p)Q^t1+pR^t1) = b/a
def LL(Q,R,l,cutoff,t1,t2):
  Q1=Q**(t1/mgt)
  Q2=Q**(t2/mgt)
  R1=R**(t1/mgt)
  R2=R**(t2/mgt)
  tot=0
  for (r,n,A,B) in l:
    if n<cutoff: continue
    def L1(p):
      return r/p-((n-r)/(1-p) if r<n else 0)-(A+B)*(R1+R2-Q1-Q2)/((1-p)*(Q1+Q2)+p*(R1+R2))+A*(R1-Q1)/((1-p)*Q1+p*R1)+B*(R2-Q2)/((1-p)*Q2+p*R2)
    p0=0.5
    while L1(p0)<0: p0/=3
    if r==n:
      p1=1
    else:
      p1=0.5
      while L1(p1)>0: p1=(2+p1)/3
    while 1:
      p=(p0+p1)/2
      if p1-p0<1e-5: break
      if L1(p)>0: p0=p
      else: p1=p
    lam=(A+B)/((1-p)*(Q1+Q2)+p*(R1+R2))
    a=lam*((1-p)*Q1+p*R1)
    b=lam*((1-p)*Q2+p*R2)
    tot+=r*log(p)+(n-r)*log(1-p)-a+A*log(a)-b+B*log(b)
  return tot

def err(xx,*aa):
  return -LL(*xx,*aa)

# Find most likely reproduction numbers for B.1.1.7 and B.1.617.2
def getrepnos(l,cutoff,t1,t2):
  xx=[1,1]
  bounds=[(0.1,10),(0.1,10)]
  res=minimize(err,xx,args=(l,cutoff,t1,t2),method="SLSQP",bounds=bounds,options={"maxiter":1000})
  if not res.success: raise RuntimeError(res.message)
  return res.x

conf=0.95
z=norm.ppf((1+conf)/2)
cutoffs=[0,20]
print("Q = estimated reproduction number of B.1.1.7")
print("R = estimated reproduction number of B.1.617.2")
print("T = increase from Q to R")
print("Date = Week-ending of second week used to calculate growth in cases")
print("Cutoff = minimum number of Sanger sequences sampled, otherwise reject LTLA")
if exclude: print("Excluding LTLAs:",exclude)
print()
print("          "+''.join("      --------- Cutoff %2d ---------"%c for c in cutoffs))
print("      Date"+''.join("         Q     R    T (%3.0f%% CI    )"%(100*conf) for c in cutoffs))
for day in range(day0,day1):
  t1=day-7-sangerday
  t2=day-sangerday
  l=[]
  for ltla in ltlas:
    (r,n)=sangnum[ltla]
    if r==0: continue
    A=sum(cases[ltla][day-day0:day-day0+7])
    B=sum(cases[ltla][day-day0+7:day-day0+14])
    l.append((r,n,A,B))
  print(daytodate(day),end='')
  for cutoff in cutoffs:
    (Q,R)=getrepnos(l,cutoff,t1,t2)
    eps=1e-2
    L0=LL(Q,R-eps,l,cutoff,t1,t2)
    L1=LL(Q,R    ,l,cutoff,t1,t2)
    L2=LL(Q,R+eps,l,cutoff,t1,t2)
    # Make confidence interval based on observed Fisher information of R (Q is unlikely to want to change that much)
    Rerr=z/sqrt((-L0+L1*2-L2)/eps**2)
    T=R/Q
    Tlow=(R-Rerr)/Q
    Thigh=(R+Rerr)/Q
    print("     %5.2f %5.2f %3.0f%% (%3.0f%% - %3.0f%%)"%(Q,R,(T-1)*100,(Tlow-1)*100,(Thigh-1)*100),end='')
  print()
