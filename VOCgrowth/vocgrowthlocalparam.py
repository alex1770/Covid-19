from stuff import *
import sys
from scipy.optimize import minimize
from random import randrange,seed
from math import log
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
mindate='2021-04-10'
#reduce=ltla2ltla
reduce=ltla2region
#reduce=ltla2country

np.set_printoptions(precision=3,linewidth=120)

lastweek=datetoday(max(sanger['WeekEndDate']))
minday=datetoday(mindate)
minday+=(lastweek-minday)%7
mindate=daytodate(minday)
nweeks=(lastweek-minday)//7

ltlas=set()
cases={}
for (date,ltla,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
  if date>mindate:
    place=reduce[ltla]
    if place not in cases: cases[place]=np.zeros([nweeks,3],dtype=int)
    week=(datetoday(date)-minday-1)//7;assert week>=0 and week<nweeks
    if var=="B.1.617.2": cases[place][week][0]+=n
    cases[place][week][1]+=n
places=sorted(list(cases))

for (ltla,date,n) in zip(apicases['areaCode'],apicases['date'],apicases['newCasesBySpecimenDate']):
  if ltla not in reduce or ltla in exclude: continue
  place=reduce[ltla]
  if place not in places: continue
  week=(datetoday(date)-minday-1)//7
  if week>=0 and week<nweeks:
    cases[place][week][2]+=n

for x in cases: print(x);print(cases[x]);print()
poi
if 1:
  for place in places:
    print("%-15s"%place,cases[place],sangnum[place])
    #for day in range(day0,day1): print(daytodate(day),cases[place][day-day0])
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
