# This is shorter/simpler than it looks. Most of the code is unused - just kept in for experimental purposes
import os,json,sys,datetime,pytz,requests,argparse
from subprocess import Popen,PIPE
from math import sqrt,log,exp
from scipy.optimize import minimize
import numpy as np
from stuff import *

np.set_printoptions(precision=4,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=10000)

infinity=7# Assume cases stabilise after this many days have elapsed
indprev=14# Assume prevalence = the sum of the last <this many> days of incidence
car=0.4#    Case ascertainment rate: assume we find this proportion of infections

parser=argparse.ArgumentParser()
parser.add_argument('-l', '--location', default='England', type=str, help='Set location: currently only England and London supported')
parser.add_argument('-s', '--skipdays', default=1, type=int, help='Discard this many days of specimen data')
parser.add_argument('-b', '--backdays', default=0, type=int, help='Do it from the point of view of this many days in the past')
args=parser.parse_args()

# Need to know public holidays (England), because they will be treated like Sundays
holidays=[datetoday(x) for x in ["2021-01-01","2021-04-02","2021-04-05","2021-05-03","2021-05-31","2021-08-30","2021-12-27","2021-12-28"]]
monday=datetoday('2021-09-20')# Any Monday
location=args.location
if location in ['England','Scotland','Wales','Northern Ireland']: areatype='nation'
else: areatype='region'

cachedir='apidata_allcaseages'
if location!='England': cachedir+='_'+location

#specmode="TimeToPublishAdjustment"
specmode="TTPadjrunningweekly"
#specmode="TTPadjdailycompound"
#specmode="Learning"
#specmode="SimpleRestrict"
#specmode="ByPublish"

#weekdayfix="NoWeekdayAdj"
#weekdayfix="SimpleAverage"
weekdayfix="MinSquareLogRatios"
#weekdayfix="MagicDeconv"

# Hospitalisation rate per 100k per day
hosprates={
  (0,10): 27.2,
  (10,20): 13.3,
  (20,30): 18.6,
  (30,40): 18.2,
  (40,50): 18.8,
  (50,60): 22.7,
  (60,70): 31.7,
  (70,80): 55.5,
  (80,90): 113.8,
  (90,150): 178.5
}
displayages=sorted(list(hosprates))
hosprates1=[hosprates[a]/1e5 for a in displayages]

def prod(l):
  p=1.
  for x in l: p*=x
  return p

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

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

d=datetime.datetime.now(pytz.timezone("Europe/London"))
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<16+10/60: today-=1# Dashboard/api updates at 4pm UK time
today-=int(args.backdays)

#minday=datetoday('2021-06-01')
minday=today-120

skipdays=int(args.skipdays)
if specmode=="ByPublish": skipdays=0

origages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
astrings=["%d_%d"%a for a in origages]

# Target save format is
# filename=publishdate, td[sex][specimendate][agerange] = cumulative cases,
# having converted agerange to open-closed format and eliminated superfluous ranges, but kept as a string because json can't handle tuples
# Note that specimendate goes back to the dawn of time, whatever minday is, because we want to save everything.
# Collect dd[publishdate]=td, td:sex -> specdate -> agestring -> number_of_cases
dd={}
os.makedirs(cachedir,exist_ok=True)
for day in range(minday-1,today+1):
  date=daytodate(day)
  fn=os.path.join(cachedir,date)
  if os.path.isfile(fn):
    with open(fn,'r') as fp: td=json.load(fp)
  else:
    male=get_data('areaType='+areatype+'&areaName='+location+'&metric=maleCases&release='+date)
    female=get_data('areaType='+areatype+'&areaName='+location+'&metric=femaleCases&release='+date)
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
          if a in origages:
            td[sexname][specdate]["%d_%d"%a]=y['value']
    with open(fn,'w') as fp: json.dump(td,fp,indent=2)
    print("Retrieved api data at",date)
  dd[date]=td

# Convert to numpy array
# ee[publishday - minday][sex 0=m, 1=f][specimenday - (minday-1)][index into displayages] = publishday's version of cumulative cases up to specimen day

reduceages={}
for (a,astr) in enumerate(astrings):
  for (da,dat) in enumerate(displayages):
    if origages[a][0]>=dat[0] and origages[a][1]<=dat[1]: reduceages[astr]=da;break
nages=len(displayages)

npub=today-minday+1
nspec=today-minday-skipdays
ee=np.zeros([npub+1,2,nspec+1,nages],dtype=int)
smindate=daytodate(minday-1)# Prepare this to compare strings because datetoday is slow
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
            if astring in reduceages:
              ee[pday][s][sday][reduceages[astring]]+=dd[pubdate][sex][specdate][astring]

# Sum over sex
# ff[publishday - (minday-1)][specimenday - (minday-1)][index into displayages] = publishday's version of cumulative cases up to specimen day
ff=ee.sum(axis=1)

# Convert from cumulative cases to new cases
# gg[publishday - (minday-1)][specimenday - minday][index into displayages] = publishday's version of new cases on specimen day
gg=ff[:,1:,:]-ff[:,:-1,:]
for i in range(nspec): gg[i+1,i,:]=0

# Convert from total up to publish day, to newly published this day
# hh[publishday - minday][specimenday - minday][index into displayages] = new cases on specimen day that were first reported on publish day
hh=gg[1:,:,:]-gg[:-1,:,:]

# Convert into (spec day, delay) co-ords
# jj[specimenday - minday][publishday-(specimenday+1)][index into displayages] = new cases on specimen day that were first reported on publish day
jj=np.zeros([npub-infinity,infinity,nages],dtype=int)
for i in range(npub-infinity):
  jj[i,:,:]=hh[i+1:i+infinity+1,i,:]

sp=np.zeros([nspec,nages],dtype=float)

# Get effective day of week (0=Monday,...,6=Sunday), where public holidays get treated as Sunday
def dow(i):
  d=minday+i
  if d in holidays: return 6
  return (d-monday)%7

# Try to undo the effect of delay from specimen to published test result
if specmode=="TimeToPublishAdjustment":
  # Use time-to-publish, based on day of week, to re-estimate full no. of cases by spec date.
  # Use more publishing days (up to "infinity") where possible; extrapolate where not (the recent days)

  # These 'time-to-publish' stats don't seem to depend that much on time of year (schools open/closed)
  # ttp[dow,dtp-1] = proportion of cases with specimen day-of-week = dow (counting from minday) , and dtp days to publish
  #                  0<=dow<7, 1<=dtp<=infinity
  ttp=np.zeros([7,infinity])
  for i in range(npub-infinity):
    d=dow(i)
    ttp[d]+=jj[i,:,:].sum(axis=1)
  for d in range(7):
    #print("%9.1f"%(t.sum()),t/t.sum()*100)
    ttp[d]/=ttp[d].sum()
  for i in range(nspec):
    n=min(npub-(i+1),infinity)
    sp[i]=hh[i+1:i+n+1,i,:].sum(axis=0)/ttp[dow(i)][:n].sum()
    #print(hh[i+1:i+n+1,i,:].sum(),ttp[dow(i)][:n].sum())
elif specmode=="TTPadjrunningweekly":
  for i in range(nspec):
    n=min(npub-(i+1),infinity)
    if n==infinity: sp[i]=hh[i+1:i+n+1,i,:].sum(axis=0)
    else: sp[i]=hh[i+1:i+n+1,i,:].sum(axis=0)/hh[i-7+1:i-7+n+1,i-7,:].sum(axis=0)*hh[i-7+1:i-7+infinity+1,i-7,:].sum(axis=0)
elif specmode=="TTPadjdailycompound":
  for i in range(nspec):
    n=min(npub-(i+1),infinity)
    if n==infinity: sp[i]=hh[i+1:i+n+1,i,:].sum(axis=0)
    else:
      sp[i]=hh[i+1:i+n+1,i,:].sum(axis=0)
      for j in range(1,infinity-n+1):
        sp[i]=sp[i]*hh[i-j+1:i+n+1,i-j,:].sum(axis=0)/hh[i-j+1:i+n,i-j,:].sum(axis=0)
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

sp[-indprev:,:].sum(axis=0)/car
for i in range(indprev,nsamp+1):
  samp=sp[i-indprev:i,:].sum(axis=0)/car
  withcov=(samp*hosprates1).sum()
  print(Date(minday+i-1)," ",location,"%5.0f"%withcov)
  
