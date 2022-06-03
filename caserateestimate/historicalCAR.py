import requests,sys
import numpy as np
from math import log
from stuff import *

mindate="2021-01-01"
if len(sys.argv)>1: mindate=sys.argv[1]
print("Mindate =",mindate)

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(10):
    try:
      response = requests.get(url+req, timeout=5)
      if response.ok: break
      error=response.text
    except BaseException as err:
      error=str(err)
  else: raise RuntimeError('Request failed: '+error)
  return response.json()['body'][::-1]


population=55.98e6
c=loadcsv("2022-05-27-CIS-weeklyprev-England-edited.csv")

now=Date(datetime.datetime.now().strftime("%Y-%m-%d"))

onsprev=[]
for (daterange,percent) in zip(c["Time period"],c["PercentPrevalence"]):
  dr=daterange.split(" to ")
  onsprev.append((Date(dr[0]),Date(dr[1]),percent/100*population))

key="newCasesBySpecimenDate"
data=get_data("areaType=nation&areaName=England&metric="+key)
date0=Date(data[0]["date"])
cases=[item[key] for item in data]
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

def pr(offset,duration):
  for item in onsprev:
    X,Y,c=item
    if X<mindate: continue
    onsest=c
    est=(cum2cases[Y-date0-offset]-cum2cases[X-date0-1-offset] - cum2cases[Y-date0-offset-duration] + cum2cases[X-date0-1-offset-duration])/(Y-X+1)
    print(X,Y,est/c)

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

pr(offset,duration)
