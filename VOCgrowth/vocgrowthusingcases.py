from stuff import *
import sys
from scipy.optimize import minimize
from random import randrange,seed
from math import log
import numpy as np

# Get ltla_2021-05-17.csv from https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv
# Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv

apicases=loadcsv("ltla_2021-05-17.csv")
sanger=loadcsv("lineages_by_ltla_and_week.2021-05-08.tsv",sep='\t')
ltladata=loadcsv("Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv")
ltla2region=dict(zip(ltladata['LAD19CD'],ltladata['RGN19NM']))
ltla2ltla=dict(zip(ltladata['LAD19CD'],ltladata['LAD19CD']))
ltla2country=dict((ltla,"England") for ltla in ltladata['LAD19CD'])

mgt=5# Mean generation time in days
mindate='2021-04-15'
minvarcases=10
hist=20
removelast=3
#reduce=ltla2ltla
reduce=ltla2region
#reduce=ltla2country

np.set_printoptions(precision=3,linewidth=120)

# Possibly use the penultimate week instead, as the last week is incomplete
sangerdate=daytodate(datetoday(max(sanger['WeekEndDate']))-7)

ltlas=set()
sangnum={}
for (date,ltla,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
  if date==sangerdate:
    place=reduce[ltla]
    if place not in sangnum: sangnum[place]=np.zeros(2,dtype=int)
    if var=="B.1.617.2": sangnum[place][0]+=n
    sangnum[place][1]+=n
places=sorted(list(place for place in sangnum if sangnum[place][0]>=minvarcases))

sangerday=datetoday(sangerdate)
day0=sangerday-hist
day1=datetoday(max(apicases['date']))-removelast+1
ndays=day1-day0
cases={}
for (ltla,date,n) in zip(apicases['areaCode'],apicases['date'],apicases['newCasesBySpecimenDate']):
  if ltla not in reduce: continue
  place=reduce[ltla]
  if place not in places: continue
  day=datetoday(date)
  if day>=day0 and day<day1:
    if place not in cases: cases[place]=np.zeros(ndays,dtype=int)
    cases[place][day-day0]+=n

if 1:
  for place in places:
    print("%-15s"%place,cases[place],sangnum[place])
  print()

# Use date <-> week end convention, compatible with Sanger
t0=sangerday-day0
for place in places:
  print(place)
  Ra=0.8
  #Rb=1.4
  (nvar,ntot)=sangnum[place]
  p=nvar/ntot
  ca=cases[place]
  # Model ca[t] as ( (1-p)*Ra^((t-t0)/mgt) + p*Rb^((t-t0)/mgt) )*ca[t0]
  c0=sum(ca[t0-6:t0+1])/7
  for t in range(6,ndays):
    c=sum(ca[t-6:t+1])/7
    #e=((1-p)*Ra**((t-t0)/mgt)+p*Rb**((t-t0)/mgt))*c0
    #print(daytodate(day0+t)," %5.1f  %5.1f"%(c,e),end='')
    print(daytodate(day0+t)," %5.1f"%c,end='')
    x=(c/c0-(1-p)*Ra**((t-t0)/mgt))/p
    if x>0:
      print("  %5.2f"%log(x),end='')
    else:
      print("  -----",end='')
    if t==t0: print(" *")
    else: print()
  print()
