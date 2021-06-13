import json,os,time,calendar,requests,sys
import numpy as np
from math import sqrt
from scipy.optimize import minimize
from scipy.stats import gamma as gammadist
from scipy.special import gammaln

if len(sys.argv)==1:
  locationsize="region"
else:
  locationsize="ltla"
  
if locationsize!="ltla" and locationsize!="region": raise RuntimeError("Unknown locationsize: "+locationsize)

#mindate='2020-10-01'
mindate='2021-03-01'

print("Using locations:",locationsize)

np.set_printoptions(precision=3,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=100000)

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def apicall(params):
  endpoint="https://api.coronavirus.data.gov.uk/v2/data"
  response = requests.get(endpoint, params=params, timeout=20)
  json=response.json()
  if not response.ok: raise RuntimeError(json["messsage"])
  return json["body"]

def getdata_ltla(mindate):
  cdir="casedata"
  os.makedirs(cdir,exist_ok=True)
  laststored=max((x for x in os.listdir(cdir) if x[:2]=="20"),default="0000-00-00")
  # Small test query to get the latest date
  data=apicall({"areaType":"ltla", "areaName":"Barnet", "metric":"newCasesBySpecimenDate", "format":"json"})
  lastavail=data[0]['date']
  if lastavail>laststored:
    data=apicall({"areaType":"ltla", "metric":"newCasesBySpecimenDate", "format":"json"})
    with open(os.path.join(cdir,lastavail),'w') as fp:
      json.dump(data,fp,indent=2)
  else:
    with open(os.path.join(cdir,laststored),'r') as fp:
      data=json.load(fp)

  # Load example map file to get location data
  with open('zoemapdata/2021-01-17','r') as fp: zm=json.load(fp)
  
  zoepl=dict((pl['lad16cd'],pl['lad16nm']) for pl in zm.values())
  offpl=dict((x['areaCode'],x['areaName']) for x in data)
  
  if 0:
    print("Different names, same lad16cd")
    for pl in zoepl:
      if pl in offpl:
        if zoepl[pl]!=offpl[pl]: print(pl,zoepl[pl],"/",offpl[pl])
    print()
    
    print("lad16cd in Zoe, not in Official")
    for pl in zoepl:
      if pl not in offpl: print(pl,zoepl[pl])
    print()
    
    print("lad16cd in Official, not in Zoe")
    for pl in offpl:
      if pl not in zoepl: print(pl,offpl[pl])
    print()
  
  # For the moment only use the places that are consistent between Zoe and Official
  # alter
  locs=set(zoepl).intersection(set(offpl))
  locs.remove('E09000012')# Hackney vs Hackney and City of London

  dates=sorted(list(set(x['date'] for x in data if x['date']>=mindate)))
  lastdate=dates[-1]
  both=set((x['date'],x['areaCode']) for x in data if x['date']>=mindate)
  r=0
  for date in dates:
    for loc in list(locs):
      if (date,loc) not in both and loc in locs: locs.remove(loc);r+=1
  if r>0: print("Removed %d locations that didn't have up to date case counts"%r)
  
  locs=sorted(list(locs))
  locind=dict((loc,i) for (i,loc) in enumerate(locs))
  
  # Convert Zoe and Official data into convenient sequences
  # zvals[ location, day ]     N x n
  # cases[ location, day ]
  N=len(locs)
  minday=datetoday(mindate)
  n=datetoday(lastdate)-minday+1
  zvals=np.zeros([N,n])# (predicted cases)/(respondents)*(population)  # , (corrected covid positive)
  cases=np.zeros([N,n],dtype=int)
  zdates=np.zeros(n,dtype=bool)
  for d in data:
    t=datetoday(d['date'])-minday
    if t>=0 and d['areaCode'] in locind: cases[locind[d['areaCode']],t]=d['newCasesBySpecimenDate']
  for date in os.listdir('zoemapdata'):
    if date[:2]=='20' and date>=mindate and date<=lastdate:
      with open(os.path.join('zoemapdata',date),'r') as fp:
        zm=json.load(fp)
        t=datetoday(date)-minday
        zdates[t]=1
        for d in zm.values():
          if d['lad16cd'] in locind:
            i=locind[d['lad16cd']]
            zvals[i,t]=d['predicted_covid_positive_count']/d['respondent']*d['population']
            #zvals[i,t]=d['corrected_covid_positive']
  # Interpolate missing dates
  for t in range(n):
    if not zdates[t]:
      for t0 in range(t-1,-1,-1):
        if zdates[t0]: break
      else: raise RuntimeError("Missing Zoe entry at start")
      for t1 in range(t+1,n):
        if zdates[t1]: break
      else: raise RuntimeError("Missing Zoe entry at end")
      zvals[:,t]=((t1-t)*zvals[:,t0]+(t-t0)*zvals[:,t1])/(t1-t0)
  
  return [offpl[loc] for loc in locs],dates,zvals,cases
            

def getdata_region(mindate):
  apidata=apicall({"areaType":"region", "metric":"newCasesBySpecimenDate", "format":"json"})
  
  # Load example map file to get location data
  with open('zoemapdata/2021-01-17','r') as fp: zm=json.load(fp)
  
  zoepl=set(pl['region'] for pl in zm.values() if pl['country']=='England')
  offpl=set(x['areaName'] for x in apidata)
  assert zoepl==offpl
  locs=sorted(list(offpl))
  locind=dict((loc,i) for (i,loc) in enumerate(locs))
  dates=sorted(list(set(x['date'] for x in apidata if x['date']>=mindate)))
  lastdate=dates[-1]
  
  # Convert Zoe and Official data into convenient sequences
  # zvals[ location, day ]     N x n
  # apicases[ location, day ]
  N=len(locs)
  minday=datetoday(mindate)
  n=datetoday(lastdate)-minday+1
  zvals=np.zeros([N,n])# (predicted cases)/(respondents)*(population)  # , (corrected covid positive)
  apicases=np.zeros([N,n],dtype=int)
  zdates=np.zeros(n,dtype=bool)
  for d in apidata:
    t=datetoday(d['date'])-minday
    if t>=0 and d['areaName'] in locind: apicases[locind[d['areaName']],t]=d['newCasesBySpecimenDate']
  for date in os.listdir('zoemapdata'):
    if date[:2]=='20' and date>=mindate and date<=lastdate:
      with open(os.path.join('zoemapdata',date),'r') as fp:
        zm=json.load(fp)
        t=datetoday(date)-minday
        zdates[t]=1
        for d in zm.values():
          if d['country']=='England' and d['region'] in locind:
            i=locind[d['region']]
            zvals[i,t]+=d['predicted_covid_positive_count']/d['respondent']*d['population']
            #zvals[i,t]+=d['corrected_covid_positive']
  # Interpolate missing dates
  for t in range(n):
    if not zdates[t]:
      for t0 in range(t-1,-1,-1):
        if zdates[t0]: break
      else: raise RuntimeError("Missing Zoe entry at start")
      for t1 in range(t+1,n):
        if zdates[t1]: break
      else: raise RuntimeError("Missing Zoe entry at end")
      zvals[:,t]=((t1-t)*zvals[:,t0]+(t-t0)*zvals[:,t1])/(t1-t0)
  
  return locs,dates,zvals,apicases

print("Loading data")
if locationsize=="ltla": locs,dates,zvals,apicases=getdata_ltla(mindate)
elif locationsize=="region": locs,dates,zvals,apicases=getdata_region(mindate)
else: assert 0  
print("Done")


N=len(locs)
n=len(dates)
assert zvals.shape==apicases.shape==(N,n)

for offset in range(30):
  a=np.array(zvals[:,:n-offset],copy=True).reshape((-1,))
  b=np.array(apicases[:,offset:],copy=True,dtype=float).reshape((-1,))
  a-=a.sum()/(N*(n-offset))
  b-=b.sum()/(N*(n-offset))
  r=np.dot(a,b)/sqrt(np.dot(a,a)*np.dot(b,b))
  #print(offset,r)

# Can see the scale factor changes. First ~70 days totcases ~= totsymp/16. Last ~100 days totcases ~= totsymp/10.
#for t in range(n):
#  print(zvals[:,t].sum(),apicases[:,t].sum())

import PtoI_noncheat

if 1:
  region='England'
  a=zvals.sum(axis=0);b=apicases.sum(axis=0)
else:
  region='North West'
  i=locs.index(region);a=zvals[i];b=apicases[i]

kernel=PtoI_noncheat.getkernel()
ad=PtoI_noncheat.deconvolve_noncheat(a,kernel,50)
#ad=PtoI_noncheat.deconvolve(a,kernel,5)

with open('temp','w') as fp:
  s=sum(b)
  c=[sum(b[i::7])/s*7 for i in range(7)]
  for t in range(n):
    print(dates[t],ad[t],a[t],b[t],sum(b[max(t-3,0):t+4])/len(b[max(t-3,0):t+4]),sum(b[max(t-6,0):t+1])/len(b[max(t-6,0):t+1]),b[t]/c[t%7],file=fp)

print("""gnuplot
set timefmt "%Y-%m-%d"
set format x "%Y-%m-%d"
set xdata time""")
print('plot "temp" u 1:2 w lines lw 3 title "%s: Zoe symptoms, non-cheat", "temp" u 1:6 w lines lw 3 title "%s: cases, 7 day lagged average"'%(region,region))

# ad[] is about 7 days ahead of sum(b[t-7:t])/7
# but that smoothing of b loses 3.5 days, and there might be a clever smoothing that doesn't lose any.
# But that would really be pushing it, and in any case at least 3.5 days ahead.
# I wonder if Zoe symptom data also has a weekday effect that could be corrected for. Doesn't look like it.

"""
if 0:
  # Partial-cheat rescaling by time-independent location-dependent function
  r=zvals.sum(1)/(cases[:,:].sum(1))
  zvals=zvals/r[:,np.newaxis]
  # Helps ccp-based estimate. Doesn't seem to help pcpc-based estimate.

def LL(zvals,cases,kern,ahead):
  k=kern.shape[-1]
  ll=0
  if deb: print("Using kern",kern)
  for l in range(N):
    lam=np.convolve(zvals[l,:n-ahead],kern[::-1],'valid')+1e-2
    targ=cases[l][ahead+k-1:]
    # lam predicts targ
    # ll+=-np.inner(targ-lam,targ-lam)
    ll+=(-lam+targ*np.log(lam)-gammaln(targ+1)).sum()
    if deb: print(l,ll)
    if deb>=2: print(lam);print(targ);print()
  return ll

def LLgamma(zvals,cases,kernsize,ahead,shape,scale):
  l=np.append(gammadist.cdf(range(kernsize),shape,scale=scale),1)
  kern=l[1:]-l[:-1]
  return LL(zvals,cases,kern,ahead)

def err(xx):
  shape,scale=xx
  return -LLgamma(zvals,cases,kernsize,ahead,shape,scale)/(N*(n-ahead-kernsize+1))

deb=0
zvals=zvals*(cases.sum()/zvals.sum())
kernsize=5
for ahead in range(20):
  res=minimize(err,[3,3],method="SLSQP",bounds=[(0.1,10),(1,20)])#,options={'ftol':1e-9,'eps':1e-6})
  if not res.success: raise RuntimeError(res.message)
  print(kernsize,ahead,res.fun,res.x)
  
# Not working well at the moment
# Can have lam=0 trying to predict case>0
# Neg bin would not fix this.
# Could either do over a larger scale (and use neg bin), OR
# Do MCMC with latent actual incidence vector - though this has zillions more parameters.

"""
