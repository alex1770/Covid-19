import os,json,sys,datetime,pytz
from requests import get
from subprocess import Popen,PIPE
from math import sqrt,log,exp
from scipy.optimize import minimize
import numpy as np
from stuff import *

np.set_printoptions(precision=4,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=10000)

cachedir='apidata_allcaseagesregion'
infinity=7# Assume cases stabilise after this many days have elapsed

# Need to know public holidays (England), because they will be treated like Sundays
holidays=[datetoday(x) for x in ["2021-01-01","2021-04-02","2021-04-05","2021-05-03","2021-05-31","2021-08-30","2021-12-27","2021-12-28"]]
monday=datetoday('2021-09-20')# Any Monday

#specmode="TimeToPublishAdjustment"
specmode="TTPadjrunningweekly"
#specmode="SimpleRestrict"
#specmode="ByPublish"

#weekdayfix="SimpleAverage"
weekdayfix="MinSquareLogRatios"
#weekdayfix="MagicDeconv"

# These need to be disjoint at the moment
displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,65),(65,150)]
#displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,40),(40,50),(50,150)]
#displayages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
#displayages=[(a,a+10) for a in range(0,80,10)]+[(80,150)]
#displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,35),(35,50),(50,65),(65,150)]
#displayages=[(0,5),(5,10),(10,15),(15,20),(20,25),(25,65),(65,80),(80,150)]

regions=['East Midlands', 'East of England', 'London', 'North East', 'North West', 'South East', 'South West', 'West Midlands', 'Yorkshire and The Humber']

# ONS 2020 population estimates from https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/09/COVID-19-weekly-announced-vaccinations-16-September-2021.xlsx
ONSpop={
  (0,18): 12093288,
  (18,25): 4709589,
  (25,30): 3771493,
  (30,35): 3824652,
  (35,40): 3738209,
  (40,45): 3476303,
  (45,50): 3638639,
  (50,55): 3875351,
  (55,60): 3761782,
  (60,65): 3196813,
  (65,70): 2784300,
  (70,75): 2814128,
  (75,80): 2009992,
  (80,150): 2855599,
}

# Simple interpolation of lower ages (which are fairly uniform per year)
ONSpop[(0,5)]=ONSpop[(0,18)]*5/18
ONSpop[(5,10)]=ONSpop[(0,18)]*5/18
ONSpop[(10,15)]=ONSpop[(0,18)]*5/18
ONSpop[(15,20)]=ONSpop[(0,18)]*3/18+ONSpop[(18,25)]*2/7
ONSpop[(20,25)]=ONSpop[(18,25)]*5/7
# Estimate older ages from typical relative proportions in 80-85, 85-90, 90+
ONSpop[(80,85)]=ONSpop[(80,150)]*254/(254+158+95)
ONSpop[(85,90)]=ONSpop[(80,150)]*158/(254+158+95)
ONSpop[(90,150)]=ONSpop[(80,150)]*95/(254+158+95)

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
if d.hour+d.minute/60<16+15/60: today-=1# Dashboard/api updates at 4pm UK time

#minday=datetoday('2021-06-01')
minday=today-60

skipdays=1
if specmode=="ByPublish": skipdays=0

origages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
astrings=["%d_%d"%a for a in origages]

# Target save format is
# filename=publishdate, td[region][sex][specimendate][agerange] = cumulative cases,
# having converted agerange to open-closed format and eliminated superfluous ranges, but kept as a string because json can't handle tuples
# Note that specimendate goes back to the dawn of time, whatever minday is, because we want to save everything.
# Collect dd[publishdate]=td, td:region -> sex -> specdate -> agestring -> number_of_cases
dd={}
os.makedirs(cachedir,exist_ok=True)
for day in range(minday-1,today+1):
  date=daytodate(day)
  fn=os.path.join(cachedir,date)
  if os.path.isfile(fn):
    with open(fn,'r') as fp: td=json.load(fp)
  else:
    male=get_data('areaType=region&metric=maleCases&release='+date)
    female=get_data('areaType=region&metric=femaleCases&release='+date)
    td={reg:{} for reg in regions}
    for sex in [male,female]:
      sexname=sex[0]['metric'][:-5]
      for reg in regions: td[reg][sexname]={}
      for d in sex:
        reg=d['areaName']
        specdate=d['date']
        td[reg][sexname][specdate]={}
        x=d[d['metric']]
        for y in x:
          a=parseage(y['age'])
          if a in origages:
            td[reg][sexname][specdate]["%d_%d"%a]=y['value']
    with open(fn,'w') as fp: json.dump(td,fp,indent=2)
    print("Retrieved api data at",date)
  dd[date]=td

def sanitise(reg):
  return reg.replace(' ','_')
  
for reg in regions:
    
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
    for sex in dd[pubdate][reg]:
      s=['male','female'].index(sex)
      for specdate in dd[pubdate][reg][sex]:
        if specdate>=smindate:
          sday=datetoday(specdate)-(minday-1)
          assert sday>=0
          if sday<nspec+1:
            for astring in dd[pubdate][reg][sex][specdate]:
              if astring in reduceages:
                ee[pday][s][sday][reduceages[astring]]+=dd[pubdate][reg][sex][specdate][astring]
  
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
  
  if weekdayfix=="SimpleAverage":
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
  
  title='Log_2 new confirmed cases per 100k per day in '+reg+' by age range.\\nDescription: http://sonorouschocolate.com/covid19/index.php?title=CasesByAge\\nData source: https://coronavirus.data.gov.uk/ at '+daytodate(today)+'; last specimen day: '+daytodate(today-1-skipdays)
  data=[]
  n=sm.shape[0]
  for (a,ar) in enumerate(displayages):
    sa=smoothpoisson(sm[:,a],0.03)/ONSpop_reduced[a]*1e5
    data.append({
      'title': ("%d - %d years"%(ar[0],ar[1]-1) if ar[1]<150 else "%d+ years"%ar[0]),
      'values': [(daytodate(minday+i),log(sa[i])/log(2)) for i in range(n)]
    })
  makegraph(title=title, data=data, mindate=daytodate(minday), ylabel='log_2 new cases per 100k per day', outfn='logcasesbyage_'+sanitise(reg)+'.png', extra=["set ytics 1","set key top left"])
  
  # Todo:
  # Validate parameters by seeing how well they predict "groundtruth" values (after several more days)
  
# montage * -geometry 960x640 ../fred.png
