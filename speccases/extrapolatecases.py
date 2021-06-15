from stuff import *
import json,os
from requests import get
import numpy as np
from random import random

np.set_printoptions(precision=4,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=100000)

adir='apiarchive'
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
mode=0

for dow in range(7):
  # dow: 0=Monday, ..., 6=Sunday
  for pred in [2]:#range(1,8):
    # Predicting cases from specimen day r-pred from information available on day r.
    numrel=3# Number of releases to go back
    numpd=3# Number of days to go back within release
    day0=minday+numrel+(monday+dow-minday-numrel)%7
    # Assemble inputs, output
    l=[];m=[]
    day=day0
    while day<maxday-infinity:
      # Predicting specimen day: 'day'-pred
      if mode==0:
        l.append([cases[day-minday-r][pred-1+p] for r in range(numrel) for p in range(numpd)])
        for r in range(numrel):
          for p in range(numpd):
            #print("Using spec day %s published day %s to help with spec day %s"%(daytodate(day-r-(pred+p)),daytodate(day-r),daytodate(day-pred)))
            pass
      elif mode==1:
        l.append([cases[day-minday-r][pred-1-r] for r in range(pred)])
        #l[-1].append(random())
        if 0:
          for r in range(pred):
            print("Using spec day %s published day %s to help with spec day %s"%(daytodate(day-pred),daytodate(day-r),daytodate(day-pred)))
      elif mode==2:
        l.append([cases[day-minday-r][pred-1-r+p] for r in range(numrel) for p in range(numpd)])
        for r in range(numrel):
          for p in range(numpd):
            #print("Using spec day %s published day %s to help with spec day %s"%(daytodate(day-(pred+p)),daytodate(day-r),daytodate(day-pred)))
            pass
      else: assert 0
      m.append(cases[day-minday+infinity-(pred-1)][infinity])
      if 0:
        print("GTR %s %d from spec date %s on published date %s"%(
          daytodate(day-pred),m[-1],
          daytodate(day-pred),
          daytodate(day+infinity-(pred-1))))
        print()
      day+=7
    A=np.array(l)
    B=np.array(m)
    C=np.linalg.lstsq(A,B)
    r=np.max(np.matmul(A,C[0])-B,0)
    print(C[0],C[1]/len(B),np.dot(r,r)/len(B))
  #print()
    
