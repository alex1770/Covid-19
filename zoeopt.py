import json,os,time,calendar,requests
import numpy as np
from scipy.optimize import minimize
from scipy.stats import gamma as gammadist
from scipy.special import gammaln
  
np.set_printoptions(precision=3,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=100000)

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

# https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesByPublishDate&format=json
# https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&areaName=Barnet&metric=newCasesByPublishDate&format=json

def apicall(params):
  endpoint="https://api.coronavirus.data.gov.uk/v2/data"
  response = requests.get(endpoint, params=params, timeout=20)
  json=response.json()
  if not response.ok: raise RuntimeError(json["messsage"])
  return json["body"]

def getcasedata():
  cdir="casedata"
  os.makedirs(cdir,exist_ok=True)
  laststored=max((x for x in os.listdir(cdir) if x[:2]=="20"),default="0000-00-00")
  # Small test query to get the latest date
  data=apicall({"areaType":"ltla", "areaName":"Barnet", "metric":"newCasesByPublishDate", "format":"json"})
  lastavail=data[0]['date']
  if lastavail>laststored:
    data=apicall({"areaType":"ltla", "metric":"newCasesByPublishDate", "format":"json"})
    with open(os.path.join(cdir,lastavail),'w') as fp:
      json.dump(data,fp,indent=2)
  else:
    with open(os.path.join(cdir,laststored),'r') as fp:
      data=json.load(fp)
  return data

print("Loading data")
data=getcasedata()

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
locs=set(zoepl).intersection(set(offpl))
locs.remove('E09000012')# Hackney vs Hackney and City of London

mindate='2020-10-01'

dates=set(x['date'] for x in data if x['date']>=mindate)
lastdate=max(dates)
both=set((x['date'],x['areaCode']) for x in data if x['date']>=mindate)
r=0
for date in sorted(list(dates)):
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
  if t>=0 and d['areaCode'] in locind: cases[locind[d['areaCode']],t]=d['newCasesByPublishDate']
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

print("Optimising")

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

