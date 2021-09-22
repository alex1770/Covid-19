import time,calendar,os,json,sys,datetime,pytz
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

# Treat public holidays (England) like Sundays
holidays=["2021-01-01","2021-04-02","2021-04-05","2021-05-03","2021-05-31","2021-08-30","2021-12-27","2021-12-28"]
monday=datetoday('2021-09-20')# Any Monday

specmode="TTP"
#specmode="SimpleRestrict"
#specmode="ByPublish"

smoothmode="SimpleAverage"
smoothmode="MinSquareLogRatios"

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(3):
    response = requests.get(url+req, timeout=30)
    if response.ok: break
  else: raise RuntimeError('Request failed: '+response.text)
  return response.json()['body'][::-1]

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

minday=datetoday('2021-06-01')
skipdays=0
if specmode=="ByPublish": skipdays=0

d=datetime.datetime.now(pytz.timezone("Europe/London"))
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<16+15/60: today-=1# Dashboard/api updates at 4pm UK time

ages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
astrings=["%d_%d"%a for a in ages]
nages=len(ages)

# Target save format is
# filename=publishdate, td[sex][specimendate][agerange] = cumulative cases,
# having converted agerange to open-closed format and eliminated superfluous ranges, but kept as a string because json can't handle tuples
# Note that specimendate goes back to the dawn of time, whatever minday is, because we want to save everything.
# Collect dd[publishdate]=td
dd={}
os.makedirs(cachedir,exist_ok=True)
for day in range(minday-1,today+1):
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
# ee[publishday - minday][sex 0=m, 1=f][specimenday - (minday-1)][index into ages] = publishday's version of cumulative cases up to specimen day

npub=today-minday+1
nspec=today-minday-skipdays
ee=np.zeros([npub+1,2,nspec+1,nages],dtype=int)
smindate=daytodate(minday-1)# Compare strings because datetoday is slow
for pubdate in dd:
  pday=datetoday(pubdate)-(minday-1)
  assert pday>=0
  for sex in dd[pubdate]:
    s=['male','female'].index(sex)
    for specdate in dd[pubdate][sex]:
      if specdate>=smindate:
        sday=datetoday(specdate)-(minday-1)
        assert sday>=0
        if sday<nspec+1:
          for astring in dd[pubdate][sex][specdate]:
            a=astrings.index(astring)
            ee[pday][s][sday][a]=dd[pubdate][sex][specdate][astring]

# Sum over sex
# ff[publishday - (minday-1)][specimenday - (minday-1)][index into ages] = publishday's version of cumulative cases up to specimen day
ff=ee.sum(axis=1)

# Convert from cumulative cases to new cases
# gg[publishday - (minday-1)][specimenday - minday][index into ages] = publishday's version of new cases on specimen day
gg=ff[:,1:,:]-ff[:,:-1,:]
for i in range(nspec): gg[i+1,i,:]=0

# Convert from total up to publish day, to newly published this day
# hh[publishday - minday][specimenday - minday][index into ages] = new cases on specimen day that were first reported on publish day
hh=gg[1:,:,:]-gg[:-1,:,:]
#ii=hh[:,:,-6:].sum(axis=2)

# Convert into (spec day, delay) co-ords
# jj[specimenday - minday][publishday-(specimenday+1)][index into ages] = new cases on specimen day that were first reported on publish day
jj=np.zeros([npub-infinity,infinity,nages],dtype=int)
for i in range(npub-infinity):
  jj[i,:,:]=hh[i+1:i+infinity+1,i,:]

sp=np.zeros([nspec,nages],dtype=float)

if specmode=="TTP":
  # Use time-to-publish, based on day of week, to re-estimate full no. of cases by spec date.
  # Use more publishing days (up to "infinity") where possible; extrapolate where not (the recent days)

  # These 'time-to-publish' stats don't seem to depend that much on time of year (schools open/closed)
  # ttp[dow,dtp-1] = proportion of cases with specimen day-of-week = dow (counting from minday) , and dtp days to publish
  #                  0<=dow<7, 1<=dtp<=infinity
  ttp=np.zeros([7,infinity])
  for d in range(7):
    t=jj[d::7,:,:].sum(axis=2).mean(axis=0)
    #print("%9.1f"%(t.sum()),t/t.sum()*100)
    ttp[d]=t/t.sum()

  for i in range(nspec):
    n=min(npub-(i+1),infinity)
    sp[i]=hh[i+1:i+infinity+1,i,:].sum(axis=0)/ttp[i%7][:n].sum()
    #print(hh[i+1:i+infinity+1,i,:].sum(),ttp[i%7][:n].sum())

elif specmode=="SimpleRestrict":
  for i in range(nspec):
    sp[i]=hh[i+1:i+2+skipdays,i,:].sum(axis=0)

elif specmode=="ByPublish":
  sp=hh[infinity:,:,:].sum(axis=1)
  minday+=infinity-2# 2 = est phase lag from using publish day

else: 
  raise RuntimeError('Unrecognised specmode '+specmode)

nsamp=sp.shape[0]

# Now sp[specimenday-minday][age index] = Est no. of samples. (In "ByPublish" mode, it's a bit of a fudged specimen day.)

print(sp[:,-6:].sum(axis=1))

#sp=sp[:,-6:]#alter

sps=sp.sum(axis=1)
if smoothmode=="SimpleAverage":
  dow=np.zeros(7)
  for d in range(7):
    dow[d]=sps[d::7].sum()/len(sps[d::7])
  dow*=7/dow.sum()
  for i in range(nsamp): sp[i,:]/=dow[i%7]
elif smoothmode=="MinSquareLogRatios":
  def slr(xx):
    pass
  pass
else: 
  raise RuntimeError('Unrecognised smoothmode '+smoothmode)

print(sp[:,-6:].sum(axis=1))
