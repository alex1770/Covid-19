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
prevalence[1]=prevalence[1]*engpop
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

k=makekernel()

with open('graph','w') as fp:
  for (i,x) in enumerate(prevalence[1][:,0]):
    print(daytodate(prevalence[0]+i),x,file=fp)
