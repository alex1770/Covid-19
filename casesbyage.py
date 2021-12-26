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
#specmode="TTPadjrunningweekly"
specmode="TTPadjdailyweekly"
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

reduceages={}
for (a,astr) in enumerate(astrings):
  for (da,dat) in enumerate(displayages):
    if origages[a][0]>=dat[0] and origages[a][1]<=dat[1]: reduceages[astr]=da;break
nages=len(displayages)
ONSpop_reduced=np.zeros(nages)
for (a,astr) in enumerate(astrings):
  if astr in reduceages: ONSpop_reduced[reduceages[astr]]+=ONSpop[origages[a]]

dd={}
for day in Daterange(minday-1,today+1):
  dd[day]=getcasesbyage_raw(day,location)

# Convert to numpy array, taking difference of cumulative values to get incremental values
npub,nspec0,cc,cn,nn=convcasesbyagetonumpy(dd,minday,today,ages=displayages)
nspec=nspec0-skipdays
hh=nn[:,:nspec,:,:].sum(axis=2)# Sum over sexes

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

if 0:
  # Print day 3 adjustments
  for i in range(nspec-2):
    print(daytodate(minday+i),i%7,hh2[i+3,i]/(hh2[i+1,i]+hh2[i+2,i]))
  print()

if 0:
  # Predict from the total of last 'fr' specimen days ==> total of last 'infinity' specimen days
  # - Not the best way to do it because requires 'infinity' days tail. Maybe better to predict individual day 3 from day 1, 2 etc.
  fr=skipdays+1
  data=[]
  targ=[]
  deb=0
  # hhq[publishday - minday][specimenday - minday] = new cases on specimen day that were first reported on publish day
  def getpredictors(hhq,i):
    # Predicting (for spec day i) sum of next fr publish days divided by sum of next infinity publish days
    # If i<npub-infinity+7 then we can make predictors
    # If i<npub-infinity then we have groundtruth, i.e. log(hhq[i+1:i+fr+1,i].sum()/hhq[i+1:i+infinity+1,i].sum())
    x1=log(hhq[i+1-7:i+fr+1-7,i-7].sum()/hhq[i+1-7:i+infinity+1-7,i-7].sum())
    A=B=0
    for j in range(i-7,0,-7):
      A+=hhq[j+1:j+fr+1,j].sum()
      B+=hhq[j+1:j+infinity+1,j].sum()
    x2=log(A/B)
    for j in range(i-7,0,-1):
      A+=hhq[j+1:j+fr+1,j].sum()
      B+=hhq[j+1:j+infinity+1,j].sum()
    x3=log(A/B)
    x4=0
    for j in range(infinity-fr):
      x4+=log(hhq[i-j:i+fr,i-1-j].sum()/hhq[i-j:i+fr+1,i-1-j].sum())
    return [1,x1,x2,x3,x4]
  for i in range(20,npub-infinity):
    y=log(hh2[i+1:i+fr+1,i].sum()/hh2[i+1:i+infinity+1,i].sum())# Partial sum divided by whole
    data.append(getpredictors(hh2,i))
    targ.append(y)
  if deb:
    with open('temp','w') as fp:
      for ((dum,x1,x2,x3,x4),y) in zip(data,targ):
        print("%9.6f %9.6f %9.6f %9.6f   %9.6f"%(x1,x2,x3,x4,y),file=fp)
  n=len(data[0])
  data=np.array(data)
  targ=np.array(targ)
  rr=[]
  for s in range(1,1<<n):
    l=[i for i in range(n) if (s&1<<i)]
    res=np.linalg.lstsq(data[:,l],targ)
    rr.append((l,res))
  rr.sort(key=lambda x:len(x[0])+1e-4*x[1][1][0])
  for (l,res) in rr:
    if deb: print("%15s  %9.6f"%(l,sqrt(res[1][0]/data.shape[0])),res[0])

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
elif specmode=="TTPadjdailyweekly":
  for i in range(nspec):
    n=min(npub-(i+1),infinity)
    if n==infinity: sp[i]=hh[i+1:i+n+1,i,:].sum(axis=0)
    else:
      base=[hh[i+1-r:i+n+1-r,i-r].sum(axis=0) for r in range(max(8,infinity-n))]
      targ7=hh[i+1-7:i+infinity-7,i-7].sum(axis=0)
      f0=1+sum(hh[i+n,i-r]/base[r] for r in range(1,infinity-n))
      f1=targ7/base[7]
      sp[i]=base[0]*(0.55*f0+0.45*f1)
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

dows=np.array([dow(i) for i in range(nsamp)])

# lam sets the scale on which the smoothed function can change per timestep. E.g., lam=0.03 <-> order of 3% change per timestep.
# Not currently used
def smoothSLR(seq,lam):
  logseq=np.log(seq)
  n=len(seq)
  def nll(xx):
    ll=0
    for i in range(n):
      ll+=-(logseq[i]-xx[i])**2/2
    for i in range(n-1):
      ll+=-((xx[i]-xx[i+1])/lam)**2/2
    return -ll
  res=minimize(nll,logseq,method="SLSQP",bounds=[(-10,20)]*n,options={"maxiter":10000})
  if not res.success: raise RuntimeError(res.message)
  return np.exp(res.x)

# Investigation to see what age bands should be grouped together for the purposes of weekday correction
if 0:
  dw=np.zeros([nages,7])
  for a in range(nages):
    dowweight=np.zeros(7)
    dowcount=np.zeros(7,dtype=int)
    for i in range(nsamp):
      d=dows[i]
      dowweight[d]+=sp[i,a]
      dowcount[d]+=1
    dowweight/=dowcount
    dowweight/=prod(dowweight)**(1/7)
    dw[a]=dowweight
    print("%2d %10s %5.3f"%(a,displayages[a],dowweight[0]/dowweight[2]),dowweight)
  print()
  dws=np.zeros([nages,7])
  for i in range(7):
    dws[:,i]=smoothSLR(dw[:,i],1)
    print("GH",i)
  for a in range(nages):
    print("%2d %10s %5.3f"%(a,displayages[a],dws[a,0]/dws[a,2]),dws[a])
  poi

if weekdayfix=="NoWeekdayAdj":
  sm=sp
elif weekdayfix=="SimpleAverage":
  sps=sp.sum(axis=1)
  dowweight=np.zeros(7)
  dowcount=np.zeros(7,dtype=int)
  for i in range(nsamp):
    d=dows[i]
    dowweight[d]+=sps[i]
    dowcount[d]+=1
  dowweight/=dowcount
  dowweight/=prod(dowweight)**(1/7)
  sm=sp/dowweight[dows][:,None]
elif weekdayfix=="MinSquareLogRatios":
  def slr(xx,sps):
    yy=sps/xx[dows]
    e=0
    for i in range(nsamp-1):#range(110,nsamp-1)
      e+=log(yy[i+1]/yy[i])**2
    return e
  sm=np.zeros(sp.shape)
  for a in range(nages):
    res=minimize(slr,[1]*7,args=(sp[:,a],),method="SLSQP",bounds=[(1,1)]+[(.5,1.5)]*6,options={"maxiter":10000})
    if not res.success: raise RuntimeError(res.message)
    xx=np.copy(res.x)
    xx/=prod(xx)**(1/7)
    sm[:,a]=sp[:,a]/xx[dows]
elif weekdayfix=="MagicDeconv":
  # (This is a general model that simulates delays in getting tested, but the optimiser ends up choosing a solution that
  # is pretty close to MinSquareLogRatios, so not continuing to develop this.)
  # (infectionincidence) . (get tested matrix) = sps
  # [nsamp+pad] . [nsamp+pad x nsamp] = [nsamp]
  pad=8
  # Postfix underscore indicates padded variable
  dows_=np.array([dow(i-pad) for i in range(pad+nsamp)])
  # Make get-tested-matrix: tm[pad+i,j]=probability of getting tested on day j from testable infection on day i
  # pdrop=dropout probability
  # ptest[j]=probability of getting tested on day j
  def maketestmatrix(pdrop,ptest_):
    tm=np.zeros([nsamp+pad,nsamp])
    for i in range(nsamp+pad):
      p=1# Probability of having "survived" up to this point (i.e., still infected, not tested, not given up the idea of getting tested)
      for j_ in range(i,nsamp+pad):
        if j_>=pad: tm[i,j_-pad]=p*ptest_[j_]
        p*=1-pdrop
        p*=1-ptest_[j_]
    return tm
  def getinfect_(xx,spec):
    pdrop,ptest_=xx[0],xx[1:8][dows_]
    tm=maketestmatrix(pdrop,ptest_)
    u,s,vt=np.linalg.svd(tm)
    assert len(s)==nsamp
    ut=u.transpose()
    v=vt.transpose()
    # tm = u[:,:-pad] @ np.diag(s) @ vt
    # ut[-1],...,ut[-pad] @ tm = 0
    # Was going to parametrise x s.t. x @ tm = sps  (x=infections), but actually that only affects [0:pad+1] of infect, so can just truncate
    # Ah, I see why now. The information being pushed forward is squeezed through 1 dimension - the number of "surviving" infections at a particular date.
    infect_=((spec @ v) @ np.diag(1/s)) @ ut[:-pad,:]
    return infect_
  def err(xx):
    infect_=getinfect_(xx,sps)
    infect=infect_[pad+1:]
    e=0
    for i in range(nsamp-1):
      if infect[i]<=0: e+=10+sqrt(-infect[i])
    for i in range(nsamp-2):
      if infect[i]>0 and infect[i+1]>0: e+=log(infect[i+1]/infect[i])**2
    return e
  res=minimize(err,[0.1]+[0.5]*7,method="SLSQP",bounds=[(0.01,0.9)]+[(0.01,0.99)]*7,options={"maxiter":10000})
  if not res.success: raise RuntimeError(res.message)
  sm=np.zeros([nsamp-1,nages],dtype=float)
  for a in range(nages):
    sm[:,a]=getinfect_(res.x,sp[:,a])[pad+1:]
  minday+=1
else: 
  raise RuntimeError('Unrecognised weekdayfix '+weekdayfix)

smoothmode="PseudoPoissonandSquareLogRatios"
# lam sets the scale on which the smoothed function can change per timestep. E.g., lam=0.03 <-> order of 3% change per timestep.
def smoothpoisson(seq,lam):
  #return seq
  n=len(seq)
  def nll(xx):
    ll=0
    for i in range(n):
      ll+=-exp(xx[i])+seq[i]*xx[i]
    for i in range(n-1):
      ll-=((xx[i]-xx[i+1])/lam)**2
    return -ll
  res=minimize(nll,np.log(seq),method="SLSQP",bounds=[(-10,20)]*n,options={"maxiter":10000})
  if not res.success: raise RuntimeError(res.message)
  return np.exp(res.x)

title='Estimated new confirmed cases per 100k per day in '+location+' by age range.\\nDescription: http://sonorouschocolate.com/covid19/index.php?title=CasesByAge\\nData source: https://coronavirus.data.gov.uk/ at '+daytodate(today)+'; last specimen day: '+daytodate(today-1-skipdays)
data=[]
n=sm.shape[0]
tot0=0
tot1=np.zeros(n)
for (a,ar) in enumerate(displayages):
  sa0=ONSpop_reduced[a]
  sa1=smoothpoisson(sm[:,a],0.03)
  sa=sa1/sa0*1e5
  tot0+=sa0
  tot1+=sa1
  data.append({
    'title': ("%d - %d years"%(ar[0],ar[1]-1) if ar[1]<150 else "%d+ years"%ar[0]),
    'values': [(daytodate(minday+i),sa[i]) for i in range(n)]
  })
if len(displayages)>1:
  tot=tot1/tot0*1e5
  data.append({
    'title': "Total",
    'values': [(daytodate(minday+i),tot[i]) for i in range(n)],
    'extra': 'dashtype (20,7)'
  })
makegraph(title=title, data=data, mindate=daytodate(displayminday), ylabel='New cases per 100k per day (log scale)', outfn=args.outfile, extra=["set key top left","set logscale y 2"])
