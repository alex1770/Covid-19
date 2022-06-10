import requests,sys
import numpy as np
from math import log
from stuff import *

mindate="2022-01-01"
if len(sys.argv)>1: mindate=sys.argv[1]
print("Mindate =",mindate)

population=55.98e6
c=loadcsv("2022-05-27-CIS-weeklyprev-England-edited.csv")

onsprev=[]
for (daterange,percent) in zip(c["Time period"],c["PercentPrevalence"]):
  dr=daterange.split(" to ")
  onsprev.append((Date(dr[0]),Date(dr[1]),percent/100*population))

now=apiday()
while 1:
  data=getcases_raw(now,location="England")
  if 'Bad' not in data: break
  print("Can't get api data for %s. Backtracking to most recent usable date."%now)
  now-=1
back=0
while 1:
  back+=7
  data0=getcases_raw(now-back,location="England")
  if 'Bad' not in data0: break
print("Using api data as of",now,"and comparing to",now-back)
data0={Date(d):c for (d,c) in data0.items()}
data={Date(d):c for (d,c) in data.items()}

date0=Date(min(data))
last=max(data)
# Correct recent incomplete entries, based on previous week
for d in Daterange(last-6,last+1):
  data[d]=round(data[d]*data[d-back]/data0[d-back])
cases=list(data.values())

cumcases=np.cumsum(cases)
# cumcases[i] = cases up to and including date0+i
cum2cases=np.cumsum(cumcases)
# cum2cases[i] = cumcases up to and including date0+i

# An onsprev item runs from dates X to Y inclusive. Here Y-X+1 = 7 or 14 (latterly 7)
# Hypothesis: is that 
# CAR(around X,Y) * onsprev(a random day between X and Y)
#                                   ~= averagecaseprev(X-offset to Y-offset)                                                                                               
#                                    = 1/(Y-X+1)*sum_{t=X}^Y caseprev(t-offset)
#                                    = 1/(Y-X+1)*sum_{t=X}^Y sum_{u=t-offset-duration+1}^{t-offset} cases(u)
#                                    = 1/(Y-X+1)*sum_{t=X}^Y (cumcases(t-offset)-cumcases(t-offset-duration))
#                                    = 1/(Y-X+1)*(cum2cases(Y-offset)-cum2cases(X-1-offset) - cum2cases(Y-offset-duration) + cum2cases(X-1-offset-duration))
# Y=X ==>
# CAR(around X) * onsprev(on day X)  = cumcases(X-offset) - cumcases(X-offset-duration)

def score(offset,duration):
  sc=0
  prev=None
  n=0
  for item in onsprev:
    X,Y,c=item
    if X<mindate: continue
    if X-date0-1-offset-duration<0 or Y-date0-offset>=len(cum2cases): continue
    onsest=c
    est=(cum2cases[Y-date0-offset]-cum2cases[X-date0-1-offset] - cum2cases[Y-date0-offset-duration] + cum2cases[X-date0-1-offset-duration])/(Y-X+1)
    cur=log(est/c)
    if prev!=None: sc+=(cur-prev)**2;n+=1
    prev=cur
  return sc/n

offset=-1
duration=14
print("Original  offset, duration =",offset,duration,", score =",score(offset,duration))
while 1:
  best=(1e9,)
  for o in range(offset-2,offset+3):
    for d in range(duration,duration+1):# Fixing duration FTM
      s=score(o,d)
      if s<best[0]: best=(s,o,d)
  (s,o,d)=best
  if o==offset and d==duration: break
  offset=o;duration=d
  print("Change to offset, duration =",offset,duration,", score =",s)

carlist=[]
for item in onsprev:
  X,Y,c=item
  if X<mindate: continue
  onsest=c
  est=(cum2cases[Y-date0-offset]-cum2cases[X-date0-1-offset] - cum2cases[Y-date0-offset-duration] + cum2cases[X-date0-1-offset-duration])/(Y-X+1)
  print(X,Y,est/c)
  carlist.append(((X+Y)/2,est/c))
carlist.append((int(last)+1,est/c))

fn='estdailyinfectionsEngland'
with open(fn,'w') as fp:
  date1=Date(int(carlist[0][0]+.999))
  adjcases=weekdayadj(cases[date1-date0:])
  i=0
  for date in Daterange(date1,last+1):
    if date>=carlist[i+1][0]: i+=1
    # carlist[i][0] <= date < carlist[i+1][0]
    (d0,c0)=carlist[i]
    (d1,c1)=carlist[i+1]
    c=((date-d0)*c1+(d1-date)*c0)/(d1-d0)
    print(date,"%8.1f"%(adjcases[date-date1]/c),file=fp)
print("Wrote estimated daily infections to file:",fn)
