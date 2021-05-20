from stuff import *
import sys
from scipy.optimize import minimize
from random import randrange,seed
from math import log

# Estimate growth advantage of B.1.617.2 over B.1.1.7 by correlating, over LTLAs, the growh over
# a pair of weeks with the relative prevalence of B.1.617.2 as estimated by Sanger sequencing data.

# Get ltla.csv from https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv
# Get Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv

apicases=loadcsv("ltla.csv")
sanger=loadcsv("lineages_by_ltla_and_week.2021-05-08.tsv",sep='\t')

lastsanger=max(sanger['WeekEndDate'])
#sangerday=datetoday(lastsanger)-7# Use penultimate week, because last week is incomplete
sangerday=datetoday(lastsanger)
sangerdate=daytodate(sangerday)
date0='2021-05-04'
date1=daytodate(datetoday(max(apicases['date']))-1)# Require dates to be less than date1, so this ignores last 2 days of cases-by-specimen-day
minsangerseq=0

mgt=5# Mean generation time in days
exclude=set()# set(["E08000001"])

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
ltlas=sorted(list(ltla for ltla in sangnum if sangnum[ltla][1]>=minsangerseq))

def LL(xx,l,cutoff):
  a,b=xx
  tot=0
  for (x,y,vx,vy,n) in l:
    if n>=cutoff:
      v=vy+b*b*vx
      tot+=(y-(a+b*x))**2/(2*v)+(1/2)*log(v)
  return tot

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

def getrepnos(l,cutoff,t1,t2):
  xx=[1,1]
  bounds=[(0.1,10),(0.1,10)]
  res=minimize(err,xx,args=(l,cutoff,t1,t2),method="SLSQP",bounds=bounds,options={"maxiter":1000})
  if not res.success: raise RuntimeError(res.message)
  return res.x

cutoffs=[0,20]
print("      Date"+''.join("       Cutoff %2d"%c for c in cutoffs))
data={}
for day in range(day0,day1):
  t1=day-7-sangerday
  t2=day-sangerday
  data[day]=[]
  for ltla in ltlas:
    (r,n)=sangnum[ltla]
    if r==0: continue
    A=sum(cases[ltla][day-day0:day-day0+7])
    B=sum(cases[ltla][day-day0+7:day-day0+14])
    data[day].append((r,n,A,B))
  print(daytodate(day),end='')
  for cutoff in cutoffs:
    (Q,R)=getrepnos(data[day],cutoff,t1,t2)
    print("     %5.2f %5.2f"%(Q,R),end='')
  print()
print()

print("      Date Cutoff  LTLAs      Extra transmission")
for day in range(day0,day1):
  t1=day-7-sangerday
  t2=day-sangerday
  cutoff=0
  l0=[x for x in data[day] if x[1]>=cutoff]
  n=len(l0)
  (Q0,R0)=getrepnos(l0,cutoff,t1,t2)
  # Bootstrap confidence interval
  sl=[]
  nit=100
  for it in range(nit):
    m=[l0[randrange(n)] for i in range(n)]
    (Q,R)=getrepnos(m,cutoff,t1,t2)
    sl.append(R/Q)
  sl.sort()
  conf=0.95
  T0=R0/Q0
  (Tmin,Tmax)=sl[int((1-conf)/2*nit)],sl[int((1+conf)/2*nit)]
  print(daytodate(day),"   %3d    %3d     %3.0f%% (%.0f%% CI: %3.0f - %3.0f%%)"%(cutoff,n,(T0-1)*100,conf*100,(Tmin-1)*100,(Tmax-1)*100))
