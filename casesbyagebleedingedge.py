import time,calendar,os,json,sys,datetime
from requests import get
from subprocess import Popen,PIPE
from math import sqrt,log,exp
from scipy.optimize import minimize
import numpy as np
from stuff import *

np.set_printoptions(precision=4,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=10000)

cachedir='apidata_allcaseages'
infinity=7# Assume cases stabilise after this many days have elapsed

monday=datetoday('2021-09-20')# Any Monday

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(3):
    response = requests.get(url+req, timeout=10)
    if response.ok: break
  else: raise RuntimeError('Request failed: '+response.text)
  return response.json()['body'][::-1]

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

minday=datetoday('2021-09-01')

d=datetime.datetime.utcnow()
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<16+15/60: today-=1

ages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
astrings=["%d_%d"%a for a in ages]
nages=len(ages)

# Target save format is
# filename=publishdate, td[sex][specimendate][agerange] = cumulative cases,
# having converted agerange to open-closed format and eliminated superfluous ranges, but kept as a string because json can't handle tuples
# Note that specimendate goes back to the dawn of time, whatever minday is.
# Collect dd[publishdate]=td
dd={}
os.makedirs(cachedir,exist_ok=True)
for day in range(minday,today+1):
  date=daytodate(day)
  fn=os.path.join(cachedir,date)
  if os.path.isfile(fn):
    with open(fn,'r') as fp: td=json.load(fp)
  else:
    male=get_data('areaType=nation&areaName=England&metric=maleCases&release='+date)
    female=get_data('areaType=nation&areaName=England&metric=femaleCases&release='+date)
    td={}
    for sex in [male,female]:
      sexname=sex[0]['metric'][:-5]
      td[sexname]={}
      for d in sex:
        specdate=d['date']
        td[sexname][specdate]={}
        x=d[d['metric']]
        for y in x:
          a=parseage(y['age'])
          if a in ages:
            td[sexname][specdate]["%d_%d"%a]=y['value']
    with open(fn,'w') as fp: json.dump(td,fp,indent=2)
    print("Retrieved api data at",date)
  dd[date]=td

# Convert to numpy array
# ee[publishday - minday][sex 0=m, 1=f][specimenday - (minday-infinity-1)][index into ages] = publishday's version of cumulative cases up to specimen day

# then convert from cumulative cases to new cases
# ff[publishday - minday][sex 0=m, 1=f][specimenday - (minday-infinity)][index into ages] = publishday's version of new cases on specimen day

n=today-minday
ee=np.zeros([n+1,2,n+infinity+1,nages],dtype=int)
for pubdate in dd:
  pday=datetoday(pubdate)-minday
  for sex in dd[pubdate]:
    s=['male','female'].index(sex)
    for specdate in dd[pubdate][sex]:
      sday=datetoday(specdate)-(minday-infinity-1)
      if sday>=0:
        for astring in dd[pubdate][sex][specdate]:
          a=astrings.index(astring)
          ee[pday][s][sday][a]=dd[pubdate][sex][specdate][astring]
        
ff=ee[:,:,1:,:]-ee[:,:,:-1,:]
for i in range(n):
  ff[i,:,infinity+i]=0
  
gg=ff[1:,:,:,:]-ff[:-1,:,:,:]
hh=gg.sum(axis=(1,3))

