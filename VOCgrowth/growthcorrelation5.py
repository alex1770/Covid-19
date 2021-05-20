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
sangerdate=daytodate(datetoday(lastsanger)-7)# Use penultimate week, because last week is incomplete
date0='2021-05-07'
date1='2021-05-17'
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


def err(xx,l,cutoff):
  a,b=xx
  tot=0
  for (x,y,vx,vy,n) in l:
    if n>=cutoff:
      v=vy+b*b*vx
      tot+=(y-(a+b*x))**2/(2*v)+(1/2)*log(v)
  return tot

fp=open('graph','w')
for day in range(day0,day1):
  l=[]
  for ltla in ltlas:
    (r,n)=sangnum[ltla]
    if r==0: continue
    x=r/n
    c0,c1=sum(cases[ltla][day-day0:day-day0+7]),sum(cases[ltla][day-day0+7:day-day0+14])
    if c0==0 or c1==0: continue
    y=c1/c0
    # Calculate approx variances of x and y
    p=(r+1)/(n+2);vx=p*(1-p)/n
    vy=c0/c1**3*(1+c0+c1)
    #l.append([x,y,vx,vy,n])
    l.append([x,y,vx*0,1,n])
    if day==day1-1: print("%s %5.3f %5.3f %5.3f %5.3f"%(ltla,x,y,vx,vy),file=fp)
  print(daytodate(day),end='')
  for cutoff in [0,20]:
    xx=[1,1]
    bounds=[(0.1,10),(-10,10)]
    res=minimize(err,xx,args=(l,cutoff),method="SLSQP",bounds=bounds,options={"maxiter":1000})
    if res.success:
      (a0,b0)=res.x
    else:
      raise RuntimeError(res.message)
    R=(1+b0/a0)**(mgt/7)
    print("  %5.2f"%R,end='')
  print()
"""      
      l0=[x for x in l if x[4]>=cutoff]
      n=len(l0)
      # Bootstrap confidence interval
      sl=[]
      nit=1000
      for it in range(nit):
        m=[l0[randrange(n)] for i in range(n)]
        res=minimize(err,xx,args=(m,cutoff),method="SLSQP",bounds=bounds,options={"maxiter":1000})
        if res.success:
          (a,b)=res.x
          sl.append(b)
        else:
          raise RuntimeError(res.message)
      sl.sort()
      conf=0.98
      (bmin,bmax)=sl[int((1-conf)/2*nit)],sl[int((1+conf)/2*nit)]
      print("Cutoff %2d : %3d LTLAs : Extra transmissibility %3.0f%% (%.0f%% CI: %.0f - %.0f%%)"%(cutoff,n,b0*100,conf*100,bmin*100,bmax*100))
    
"""
