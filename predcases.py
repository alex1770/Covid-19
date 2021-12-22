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

parser=argparse.ArgumentParser()
parser.add_argument('-l', '--location', default='England', type=str, help='Set location: currently only England and London supported')
parser.add_argument('-s', '--skipdays', default=1, type=int, help='Discard this many days of specimen data')
parser.add_argument('-b', '--backdays', default=0, type=int, help='Do it from the point of view of this many days in the past')
parser.add_argument('-o', '--outfile', default='logcasesbyage.png', type=str, help='Output png filename')
args=parser.parse_args()

# Need to know public holidays (England), because they will be treated like Sundays
holidays=[datetoday(x) for x in ["2021-01-01","2021-04-02","2021-04-05","2021-05-03","2021-05-31","2021-08-30","2021-12-25","2021-12-27","2021-12-28"]]
monday=datetoday('2021-09-20')# Any Monday

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

# These need to be disjoint at the moment
#displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,65),(65,150)]
#displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,40),(40,50),(50,150)]
#displayages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
#displayages=[(a,a+10) for a in range(0,80,10)]+[(80,150)]
displayages=[(a,a+10) for a in range(0,70,10)]+[(70,150)]
#displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,35),(35,50),(50,65),(65,150)]
#displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,65),(65,80),(80,150)]
#displayages=[(0,150)]
#displayages=[(0,20),(20,30),(30,50),(50,70),(70,150)]
#displayages=[(0,20),(20,30),(30,50),(50,55),(55,60),(60,65),(65,70),(70,150)]
#displayages=[(5,15)]

loclookup={
  "E92000001":"England",
  "E12000001":"North East",
  "E12000002":"North West",
  "E12000003":"Yorkshire and The Humber",
  "E12000004":"East Midlands",
  "E12000005":"West Midlands",
  "E12000006":"East of England",
  "E12000007":"London",
  "E12000008":"South East",
  "E12000009":"South West"
}

location=args.location.replace('_',' ')
location=loclookup.get(location,location)
if location in ['England','Scotland','Wales','Northern Ireland']: areatype='nation'
else: areatype='region'

cachedir='apidata_allcaseages'
if location!='England': cachedir+='_'+location

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

ONSpop={}
for (desc,acode,sex,age,n) in csvrows('ONS-population_2021-08-05.csv',['category','areaCode','gender','age','population']):
  if desc=='AGE_SEX_5YEAR' and sex=='ALL' and acode in loclookup and loclookup.get(location,location).lower()==loclookup[acode].lower():
    ONSpop[parseage(age)]=int(n)

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

d=datetime.datetime.now(pytz.timezone("Europe/London"))
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<16+5/60: today-=1# Dashboard/api updates at 4pm UK time
today-=int(args.backdays)

#minday=datetoday('2021-06-01')
minday=today-120
displayminday=today-90

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
ONSpop_reduced=np.zeros(nages)
for (a,astr) in enumerate(astrings):
  if astr in reduceages: ONSpop_reduced[reduceages[astr]]+=ONSpop[origages[a]]

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

hh2=hh[:,:].sum(axis=2)
fr=skipdays+1
nspec0=nspec+skipdays# = npub-1
# npub-fr=nspec
def clip(x,a,b): return min(max(x,a),b)
#for dayn in range(fr+1,5):#infinity):
if 1:
  preds=[]
  targs=[]
  al=0.55
  if 0:
    for spec in range(11,npub-infinity+1):
      # Not allowed to source rows past spec+fr, since that would be looking into the future
      
      #base=np.array([hh[spec+1-i:spec+fr+1-i,spec-i].sum(axis=0) for i in range(spec+1)])
      #targ=np.array([hh[spec+fr+1-i:spec+infinity-i,spec-i].sum(axis=0) for i in range(spec+1)])
      #shorttarg=np.array([hh[spec+fr+1-i,spec-i] for i in range(spec+1)])
      #targs.append(hh[spec+fr+1:spec+infinity,spec].sum(axis=0))
      
      base=np.array([hh2[spec+1-i:spec+fr+1-i,spec-i].sum(axis=0) for i in range(spec+1)])
      targ=np.array([hh2[spec+fr+1-i:spec+infinity-i,spec-i].sum(axis=0) for i in range(spec+1)])
      shorttarg=np.array([hh2[spec+fr+1-i,spec-i] for i in range(spec+1)])
      targs.append(hh2[spec+fr+1:spec+infinity,spec].sum(axis=0))
      
      #print(base[7],targ[7],base[0],targ[0])
      #preds.append(base[0])
      #preds.append(base[0]*targ[7]/base[7])
      #preds.append(base[0]*targ[1]/base[1])
      #preds.append(base[0]*targ[7]/base[7]*clip(targ[1]/base[1]/(targ[8]/base[8]),0.8,1.2))
      #print(targ[1]/base[1]/(targ[8]/base[8]))
      #preds.append(base[0]*np.sqrt(targ[1]/base[1]*targ[7]/base[7]))
      #preds.append(base[0]*(targ[1]/base[1])**al*(targ[7]/base[7])**(1-al))
      #preds.append(base[0]*targ[7::7].sum()/base[7::7].sum())
      #preds.append(base[0]*targ[1:].sum()/base[1:].sum())
      #preds.append(base[0]*targ[1:8].sum()/base[1:8].sum())
      #preds.append(base[0]*targ[7:28:7].sum()/base[7:28:7].sum())
      preds.append(base[0]*(targ[1]+targ[7])/(base[1]+base[7]))
      #preds.append(base[0]*(targ[1:8:3].sum())/(base[1:8:3].sum()))
      #preds.append(base[0]*(al*targ[1]+(1-al)*targ[7])/((1-al)*base[1]+al*base[7]))
      #preds.append(base[0]*(al*targ[1]/base[1]+(1-al)*targ[7]/base[7]))
      #preds.append(base[0]*targ[7]/base[7]*((shorttarg[1]/base[1])/(shorttarg[7]/base[7]))**al)

  if 0:
    for spec in range(11,npub-infinity+1):
      # Not allowed to source rows past spec+fr, since that would be looking into the future
      base=np.array([hh2[spec+1-i:spec+fr+1-i,spec-i].sum(axis=0) for i in range(spec+1)])
      targ=np.array([hh2[spec+fr+1-i:spec+infinity-i,spec-i].sum(axis=0) for i in range(spec+1)])
      f=sum(hh2[spec+fr,spec-i]/hh2[spec+1-i:spec+1+fr-i,spec-i].sum(axis=0) for i in range(1,infinity-fr))
      #preds.append(base[0]*f)
      #preds.append(base[0]*targ[7]/base[7])
      preds.append(base[0]*(al*f+(1-al)*targ[7]/base[7]))
      #preds.append(base[0]*sqrt(f*targ[7]/base[7]))
      targs.append(hh2[spec+fr+1:spec+infinity,spec].sum(axis=0))
  
  if 0:
    for spec in range(11,npub-infinity+1):
      # Not allowed to source rows past spec+fr, since that would be looking into the future
      base=np.array([hh[spec+1-i:spec+fr+1-i,spec-i].sum(axis=0) for i in range(spec+1)])
      targ=np.array([hh[spec+fr+1-i:spec+infinity-i,spec-i].sum(axis=0) for i in range(spec+1)])
      f=sum(hh[spec+fr,spec-i]/base[i] for i in range(1,infinity-fr))
      #preds.append(base[0]*f)
      #preds.append(base[0]*targ[7]/base[7])
      preds.append(base[0]*(al*f+(1-al)*targ[7]/base[7]))
      #preds.append(base[0]*np.sqrt(f*targ[7]/base[7]))
      targs.append(hh[spec+fr+1:spec+infinity,spec].sum(axis=0))
  
  if 0:
    for spec in range(11,npub-infinity+1):
      # Not allowed to source rows past spec+fr, since that would be looking into the future
      base=np.array([hh[spec+1-i:spec+fr+1-i,spec-i].sum(axis=0) for i in range(spec+1)])
      targ=np.array([hh[spec+1-i:spec+infinity-i,spec-i].sum(axis=0) for i in range(spec+1)])
      f=1+sum(hh[spec+fr,spec-i]/base[i] for i in range(1,infinity-fr))
      #preds.append(base[0]*f)
      #preds.append(base[0]*targ[7]/base[7])
      #print(f,targ[7]/base[7])
      preds.append(base[0]*(al*f+(1-al)*targ[7]/base[7]))
      #preds.append(base[0]*np.sqrt(f*targ[7]/base[7]))
      targs.append(hh[spec+1:spec+infinity,spec].sum(axis=0))
  
  if 1:
    for spec in range(11,npub-infinity+1):
      # Not allowed to source rows past spec+fr, since that would be looking into the future
      base=np.array([gg[spec+fr+1-i,spec-i] for i in range(spec+1)])
      targ=np.array([gg[spec+infinity-i,spec-i] for i in range(spec+1)])
      f=1+sum((gg[spec+fr+1,spec-i]-gg[spec+fr,spec-i])/base[i] for i in range(1,infinity-fr))
      preds.append(base[0]*(al*f+(1-al)*targ[7]/base[7]))
      targs.append(gg[spec+infinity,spec])
  
  n=len(preds)
  preds=np.array(preds,dtype=float)
  targs=np.array(targs,dtype=float)
  lam=(preds*targs).sum()/(preds*preds).sum()
  #lam=1
  preds*=lam
  diff=targs-preds
  if len(diff.shape)>1: diff=diff.sum(axis=1)
  print(lam,sqrt((diff@diff)/n),abs(diff).sum()/n)
  #for spec in range(7,nspec0-dayn+1):
  #  print(spec,diff[spec-7])
  #print()
