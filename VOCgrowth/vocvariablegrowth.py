from stuff import *
import sys
from scipy.optimize import minimize
from random import randrange,seed
from math import log,exp,sqrt
import numpy as np

# Get ltla.csv from https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv
# Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv

apicases=loadcsv("ltla.csv")
sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
ltladata=loadcsv("Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv")
ltla2region=dict(zip(ltladata['LAD19CD'],ltladata['RGN19NM']))
ltla2ltla=dict(zip(ltladata['LAD19CD'],ltladata['LAD19CD']))
ltla2country=dict((ltla,"England") for ltla in ltladata['LAD19CD'])
exclude=set()#set(['E08000001'])

mgt=5# Mean generation time in days
minday=datetoday('2021-04-10')# Inclusive
#reduce=ltla2ltla
reduce=ltla2region
#reduce=ltla2country

# Roughly how much daily growth rate is allowed to change in 1 day
sig=0.0025

# Case ascertainment rate
asc=0.4

np.set_printoptions(precision=3,linewidth=120)

maxday=datetoday(max(apicases['date']))-2# Inclusive
ndays=maxday-minday+1
lastsang=datetoday(max(sanger['WeekEndDate']));assert maxday>=lastsang
nweeks=(lastsang-minday+1)//7
# Sanger week number is nweeks-1-(lastsang-day)//7

sang={}
for (date,ltla,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
  day=datetoday(date)
  week=nweeks-1-(lastsang-day)//7
  if week>=0 and week<nweeks:
    place=reduce[ltla]
    if place not in sang: sang[place]=np.zeros([nweeks,2],dtype=int)
    if var=="B.1.617.2": sang[place][week][1]+=n
    else: sang[place][week][0]+=n
places=sorted(list(sang))

# Relatively stable period (inclusive) over which to take the weekday averages
weekadjdates=[datetoday('2021-04-03'),datetoday('2021-05-14')]
weekadj=np.zeros(7)
for (date,n) in zip(apicases['date'],apicases['newCasesBySpecimenDate']):
  day=datetoday(date)
  if day>=weekadjdates[0] and day<=weekadjdates[1]: weekadj[day%7]+=n
weekadjp=weekadj*7/sum(weekadj)
  
cases={}
for (ltla,date,n) in zip(apicases['areaCode'],apicases['date'],apicases['newCasesBySpecimenDate']):
  if ltla not in reduce or ltla in exclude: continue
  day=datetoday(date)
  d=day-minday
  if d<0 or d>=ndays: continue
  place=reduce[ltla]
  if place not in places: continue
  if place not in cases: cases[place]=np.zeros(ndays)
  cases[place][d]+=n/weekadjp[day%7]
#for x in places: print(x);print(cases[x]);print()

# ndays+2 parameters to be optimised:
# 0: a0
# 1: b0
# 2: h/sig
# 3 ... 3+ndays-2 : g_0/sig, ..., g_{ndays-2}/sig

def expand(xx,sig):
  (a0,b0,h)=xx[:3]
  AA=[exp(a0)];BB=[exp(b0)]
  a=a0;b=b0
  for i in range(ndays-1):
    g=xx[3+i]
    a+=g*sig;b+=(g+h)*sig
    AA.append(exp(a))
    BB.append(exp(b))
  return AA,BB

def NLL(xx,lcases,lsang,sig,p):
  #print(xx)
  AA,BB=expand(xx,sig)
  tot=0
  for i in range(ndays):
    lam=p*(AA[i]+BB[i])
    tot+=-lam+lcases[i]*log(lam)
  for i in range(ndays-2):
    tot+=-(xx[3+i+1]-xx[3+i])**2/2
  for w in range(nweeks):
    endweek=lastsang-(nweeks-1-w)*7-minday
    A=sum(AA[endweek-6:endweek+1])
    B=sum(BB[endweek-6:endweek+1])
    tot+=lsang[w][0]*log(A/(A+B))+lsang[w][1]*log(B/(A+B))
  return -tot

for place in places:
  print(place)
  print("="*len(place))
  bounds=[(-10,20),(-10,20),(-500,500)]+[(-50,50)]*(ndays-1)
  res=minimize(NLL,[0]*(ndays+2),args=(cases[place],sang[place],sig,asc),bounds=bounds,method="SLSQP",options={"maxiter":1000})
  if not res.success: raise RuntimeError(res.message)
  xx=res.x
  #print(xx)
  AA,BB=expand(xx,sig)
  (a0,b0,h)=xx[:3]
  print("A    = estimated number of new cases of non-B.1.617.2 on this day")
  print("B    = estimated number of new cases of B.1.617.2 on this day")
  print("Pred = predicted number of cases seen this day = (ascertainment rate)*(A+B)")
  print("Seen = number of cases seen this day, with weekday adjustment")
  print("Q    = estimated reproduction rate of non-B.1.617.2 on this day")
  print("R    = estimated reproduction rate of B.1.617.2 on this day")
  print("      Date       A       B    Pred    Seen       Q     R")
  for i in range(ndays):
    print(daytodate(minday+i),"%7.0f %7.0f %7.0f %7.0f"%(AA[i],BB[i],asc*(AA[i]+BB[i]),cases[place][i]),end='')
    if i<ndays-1:
      g=xx[3+i]
      print("   %5.2f %5.2f"%(exp(g*sig*mgt),exp((g+h)*sig*mgt)))
    else:
      print()
  h0=xx[2]
  ff=[0,res.fun,0]
  eps=0.01/sig
  for i in [-1,1]:
    eps=0.01/sig
    h=h0+i*eps
    xx=[0,0,h]+[0]*(ndays-1)
    bounds[2]=(h,h)
    res=minimize(NLL,[0]*(ndays+2),args=(cases[place],sang[place],sig,asc),bounds=bounds,method="SLSQP",options={"maxiter":1000})
    if not res.success: raise RuntimeError(res.message)
    ff[i+1]=res.fun
  dh=1.96/sqrt((ff[0]-2*ff[1]+ff[2])/eps**2)
  print("Transmission advantage %.0f%% (%.0f%% - %.0f%%)"%((exp(h0*sig*mgt)-1)*100,(exp((h0-dh)*sig*mgt)-1)*100,(exp((h0+dh)*sig*mgt)-1)*100))
  print()
