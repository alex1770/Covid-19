import time,calendar,os,json,sys,datetime,requests,sys
import numpy as np
from stuff import *
from scipy.stats import gamma
from math import log,exp

# ONS infection survey data is from table 1b and 1p of https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata

mintrain=datetoday('2021-01-17')
minprint=datetoday('2020-10-01')
discardspecdays=1
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
# Note that prevalence data in 2021-01-02 - 2021-01-16 are anomalous
# Also crazy behaviour before 2020-10-01
# Also other minor(?) anomalies
    
# speccasesadjust[d][r] = chance that a specimen from day of the week d (Monday=0) is reported (by) r days later
speccasesadjust=np.array([[ 0.    , 0.3954, 0.8823, 0.9664, 0.9914, 0.9975, 0.9989, 1.    ],
                          [ 0.    , 0.3264, 0.8731, 0.9707, 0.9893, 0.9962, 1.0003, 1.    ],
                          [ 0.    , 0.3561, 0.8745, 0.9604, 0.983 , 0.989 , 0.9936, 1.    ],
                          [ 0.    , 0.3511, 0.8699, 0.9567, 0.9768, 0.9873, 0.9974, 1.    ],
                          [ 0.    , 0.3132, 0.8363, 0.942 , 0.9701, 0.988 , 0.9964, 1.    ],
                          [ 0.    , 0.2943, 0.8404, 0.9218, 0.9595, 0.984 , 0.994 , 1.    ],
                          [ 0.    , 0.4585, 0.9093, 0.961 , 0.9859, 0.9934, 0.9987, 1.    ]])
specadjmax=speccasesadjust.shape[1]# Cases are considered stable at specadjmax-1 days after specimen taken
monday=datetoday('2021-06-07')# Example of a Monday

# Get cases
cases0=get_data('areaType=nation&areaName=England&metric=newCasesBySpecimenDate')
publishedday=datetoday(cases0[-1]['date'])+1# Day when api published results
lastzero=-1
for i in range(len(cases0)):
  if cases0[i]['newCasesBySpecimenDate']==0: lastzero=i
cases0=cases0[lastzero+1:]
l=[datetoday(x['date']) for x in cases0]
assert list(range(l[0],l[-1]+1))==l# Check contiguous
a=np.array([x['newCasesBySpecimenDate'] for x in cases0])
for r in range(1,specadjmax):
  dayofweek=(l[0]+len(a)-r-monday)%7
  a[-r]/=speccasesadjust[dayofweek][r]
cases=(l[0],a[:len(a)-discardspecdays])

# Make recovery kernel: kernel[i] = P(recovery time > i+0.5 days)
def makekernel(shape=2.595, scale=4.48):
  kernel=[]
  for i in range(1000000):
    x=gamma.sf(i+0.5,shape,scale=scale)
    kernel.append(x)
    if x<1e-3: break
  return np.array(kernel)

def setup(shape,scale,logcar,offset,cases,prevalence,logprev,mintrain,pr=0):
  ker=makekernel(shape,scale)
  conv=np.convolve(cases[1],ker)
  
  logpred=np.log(conv)-logcar
  # ker[i] = probability of still being +ve after i+.5 days
  # cases[1][j] = number of new confirmed cases on day cases[0]+j
  # conv[k] = number of confirmed cases that are still positive on day cases[0]+k
  # conv[k] is non-truncated if len(ker)-1 <= k < len(cases[1])
  
  delta=prevalence[0]-cases[0]+offset
  # Compare conv nominally at day d+offset with prevalence nominally at day d.
  # That means logpred[k] is compared with logprev[k+cases[0]-prevalence[0]-offset] = logprev[k-delta]
  # mintrain <= k-delta+prevalence[0] = k+cases[0]-offset
  # len(ker)-1 <= k < len(cases[1])
  # 0 <= k-delta < len(logprev)
  # ==>
  # max(mintrain-cases[0]+offset, len(ker)-1, delta) <= k < min(len(cases[1]), len(logprev)+delta)
  k0=max(mintrain-cases[0]+offset, len(ker)-1, delta)
  k1=min(len(cases[1]), len(logprev)+delta)
  return (k0,k1,delta,logpred)

def err(xx,offset,cases,prevalence,logprev,logprevsd,mintrain):
  shape,scale,logcar=xx
  (k0,k1,delta,logpred)=setup(shape,scale,logcar,offset,cases,prevalence,logprev,mintrain)
  #print(mintrain-cases[0]+offset, len(ker)-1, delta)
  #print(len(cases[1]), len(logprev)+delta)
  #print(k0,k1)
  #e=(logpred[k0:k1]-logprev[k0-delta:k1-delta])/logprevsd[k0-delta:k1-delta]
  #e=(logpred[k0:k1]-logprev[k0-delta:k1-delta])
  e=(np.exp(logpred[k0:k1])-np.exp(logprev[k0-delta:k1-delta]))/3e5
  e2=abs(xx[0]*xx[1]-11)/10
  return sum(e*e)+e2

from scipy.optimize import minimize
best=(1e100,)
for offset in range(-10,10):
  xx=[2.6, 4.5, -1]
  bounds=[(0.1,30), (0.1,30), (-5,1)]
  #bounds=[(2.6,2.6), (4.5,4.5), (-5,1)]
  res=minimize(err,xx,args=(offset,cases,prevalence,logprev,logprevsd,mintrain),bounds=bounds,method="SLSQP",options={"maxiter":10000,"eps":1e-4})
  if not res.success:
    raise RuntimeError(res.message)
  xx=res.x
  shape,scale,logcar=xx
  (k0,k1,delta,logpred)=setup(shape,scale,logcar,offset,cases,prevalence,logprev,mintrain,1)
  for x in xx:
    print('%7.3f'%x,end=' ')
  print('%3d %6.3f %7.3f %3d'%(offset,xx[0]*xx[1]+offset,res.fun,k1-k0))
  if res.fun<best[0]: best=(res.fun, offset, res.x)

val,offset,xx=best
shape,scale,logcar=xx
(k0,k1,delta,logpred)=setup(shape,scale,logcar,offset,cases,prevalence,logprev,mintrain)
ker=makekernel(shape,scale)
k2=max(minprint-cases[0]+offset, len(ker)-1, delta)
k3=max(len(cases[1]), len(logprev)+delta)

with open('graph','w') as fp:
  for k in range(k2,k3):
    l=k-delta
    print(daytodate(prevalence[0]+l),end=' ',file=fp)
    if k>=0 and k<len(cases[1]):
      print(logpred[k],end=' ',file=fp)
    else:
      print('-',end=' ',file=fp)
    if l>=0 and l<len(logprev):
      print(logprev[l],file=fp)
    else:
      print('-',file=fp)

print('Written to file "graph"')
print('gnuplot')
print('set xdata time;set timefmt "%Y-%m-%d";set format x "%Y-%m-%d"')
print('plot "graph" u 1:2 w lines lw 2 title "Cases convolved with Gamma(shape=%.3f, scale=%.3f, offset=%d) divided by CAR=%.3f", "graph" u 1:3 w lines lw 2 title "ONS prev"'%(shape,scale,offset,exp(logcar)))
