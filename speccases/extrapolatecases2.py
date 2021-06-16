from stuff import *
import json,os
from requests import get
import numpy as np
from random import random
from math import sqrt,log

np.set_printoptions(precision=4,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=100000)

adir='casesbyspecdate'
mindate='2021-02-24'# Can't get archive data earlier than 2021-02-24
now=datetime.datetime.utcnow().strftime('%Y-%m-%d')

infinity=7# Assume cases have stabilised after this many days

minday=datetoday(mindate)
maxday=datetoday(now)# exclusive
monday=datetoday('2021-06-07')

cases=[]
for day in range(minday,maxday):
  date=daytodate(day)
  fn=os.path.join(adir,date)
  if not os.path.isfile(fn):
    print("Loading cases as at",date)
    url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDate&format=json&release='+date
    response=get(url,timeout=10)
    if not response.ok: raise RuntimeError(f'Request failed: '+response.text)
    data=response.json()['body']
    with open(fn,'w') as fp: json.dump(data,fp,indent=2)
  with open(fn,'r') as fp:
    a=json.load(fp)
    l=[d['newCasesBySpecimenDate'] for d in a if d['date']!=date]# 2021-03-15 has erroneous entry for 2021-03-15
    cases.append(l)

# cases[x][y] = cases from specimen day minday+x-(y+1), as reported on day minday+x
# Specimen day minday+s, report day minday+r --> cases[r][r-s-1]

# transfer[d][r'] = number of cases on weekday d that get reported by (specimen date)+r'. r'=1, ..., infinity
transfer=np.zeros([7,infinity+1],dtype=int)
p0=np.zeros([7],dtype=int)
p1=np.zeros([7,infinity+1])
p2=np.zeros([7,infinity+1])
l1=np.zeros([7,infinity+1])
l2=np.zeros([7,infinity+1])


for s in range(len(cases)-infinity):
  d=(minday+s-monday)%7
  p0[d]+=1
  for r in range(s+1,s+infinity+1):
    transfer[d][r-s]+=cases[r][r-s-1]
    p=cases[r][r-s-1]/cases[s+infinity][infinity-1]
    p1[d][r-s]+=p
    p2[d][r-s]+=p*p
    l1[d][r-s]+=log(p)
    l2[d][r-s]+=log(p)**2

print("Across = number of days after specimen day that result is reported")
print("Down = day of the week of the specimen day, starting at Monday")
print()
print(transfer)
print()
mu0=transfer/transfer[:,infinity][:,None]
print(mu0)
print()

print("mean(p)")
mu=p1/p0[:,None]
print(mu)
print("mean(p) Python format")
print(np.array2string(mu,separator=','))

sd=np.sqrt((p2-p1**2/p0[:,None])/(p0[:,None]-1))
print()
print("sd(p)")
print(sd)
print()

#sdq=np.sqrt(mu*(1-mu)/p0[:,None])
#print(sdq)
#print()
#print(sd/sdq)

sdl=np.sqrt((l2-l1**2/p0[:,None])/(p0[:,None]-1))
print("sd(logp)")
print(sdl)
print()

