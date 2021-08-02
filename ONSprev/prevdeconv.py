import time,calendar,os,json,sys,datetime,requests,sys
import numpy as np
from stuff import *
from scipy.stats import gamma

# ONS infection survey data is from table 1b and 1p of https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata

minprev=datetoday('2021-01-17')
maxprev=datetoday('2021-07-17')
discardspecdays=2
engpop=54.53e6# This population is used by ONS as the denominator for England (age>=2)


def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(3):
    response = requests.get(url+req, timeout=10)
    if response.ok: break
  if not response.ok:
    raise RuntimeError('Request failed: '+response.text)
  return response.json()['body'][::-1]

# Need custom parser because sheet 1p is a mess
def parseONScsv():
  d2v={}
  with open('1b.csv') as fp:
    reader=csv.reader(fp)
    phase=0
    for y in reader:
      if y[0]=='Date': phase=1
      if phase==1 and y[0][:1].isdigit(): phase=2
      if phase==2 and not y[0][:1].isdigit(): phase=0
      if phase==2:
        day=datetoday(y[0])
        # Need to use number of people to get precision, because they truncate percentages to 2dp on this sheet
        d2v[day]=[float(z.replace(',',''))/engpop for z in y[4:7]]
  with open('1p.csv') as fp:
    reader=csv.reader(fp)
    phase=0
    for y in reader:
      if y[0]=='Date': phase=1
      if phase==1 and y[0][:1].isdigit(): phase=2
      if phase==2 and not y[0][:1].isdigit(): phase=0
      if phase==2:
        day=datetoday(y[0])
        # Publication dates are in reverse order, so earliest in file (latest in time) takes priority
        if day not in d2v:
          d2v[day]=[float(z[:-1])/100 for z in y[1:4]]
  l=sorted(list(d2v))
  #for day in l: print(daytodate(day),d2v[day])
  # Check contiguous
  assert list(range(l[0],l[-1]+1))==l
  return (l[0],np.array([d2v[x] for x in l]))

# Get ONS prevalence
prevalence=parseONScsv()
# Fix bug in ONS spreadsheet on 2021-05-02:
fixind=datetoday('2021-05-02')-prevalence[0]
prevalence[1][fixind]=(prevalence[1][fixind-1]+prevalence[1][fixind+1])/2
logprev=np.log(prevalence[1][:,0]*engpop)
logprevsd=np.log(prevalence[1][:,2]/prevalence[1][:,1])
# Note that data in 2021-01-02 - 2021-01-16 are anomalous
# Also crazy behaviour before 2020-10-01
# Also other minor anomalies
    
# Get cases
cases0=get_data('areaType=nation&areaName=England&metric=newCasesBySpecimenDate')
l=[datetoday(x['date']) for x in cases0]
# Check contiguous
assert list(range(l[0],l[-1]+1))==l
cases=[l[0],np.array([x['newCasesBySpecimenDate'] for x in cases0])]

# Make recovery kernel: kernel[i] = P(recovery time > i+0.5 days)
def makekernel(shape=2.595, scale=4.48):
  kernel=[]
  for i in range(1000000):
    x=gamma.sf(i+0.5,shape,scale=scale)
    kernel.append(x)
    if x<1e-3: break
  return np.array(kernel)

ker=makekernel()
con=np.convolve(cases[1],ker)
logcon=np.log(con)
# ker[i] = probability of still being +ve after i+.5 days
# cases[1][j] = number of new confirmed cases on day cases[0]+j
# con[k] = number of confirmed cases that are still positive on day cases[0]+k
# con[k] is non-truncated if len(ker)-1 <= k < len(cases[1])

offset=5
delta=prevalence[0]-cases[0]-offset
# Compare conv nominally at day d with prevalence nominally at day d+offset. (Expect offset>0.)
# That means logcon[k] is compared with logprev[k+cases[0]-prevalence[0]+offset] = logprev[k-delta]
# len(ker)-1 <= k < len(logcon)
# 0 <= k-delta < len(logprev)
# ==>
# max(len(ker)-1, delta) <= k < min(len(logcon), len(logprev)+delta)
k0=max(len(ker)-1, delta)
k1=min(len(logcon), len(logprev)+delta)

with open('graph','w') as fp:
  for k in range(k0,k1):
    l=k-delta
    print(daytodate(cases[0]+k),logcon[k],daytodate(prevalence[0]+l),logprev[l],file=fp)
