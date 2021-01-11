import csv,sys
from math import log,exp

if sys.version_info[0]<3: print("Requires Python 3");sys.exit(1)

def getcountrydata(country,thr=0,smoothlevel=0,region="",source="worldometer"):
  """
  Returns (dates, totalconfirmed, totaldeaths, recovered (if available), active (if available), newcases, newdeaths)
  for the specified country.
  Optionally set thr to ignore initial values with fewer than thr confirmed cases.
  """
  
  equivnames={}
  with open("countrynames") as fp:
    r=csv.reader(fp)
    for x in r:
      if len(x)==0 or x[0][:1]=='#': continue
      for y in x[1:]: equivnames[y]=x[0]

  country=equivnames.get(country,country)
  
  dates=[]
  confirmed=[]
  deaths=[]
  recovered=[]
  active=[]
  newc=[]
  newd=[]
  with open(source+".csv") as f:
    r=csv.reader(f)
    first=1
    for x in r:
      if first: first=0;continue
      if equivnames.get(x[1],x[1])==country:
        dates.append(x[0])
        if x[4].isdigit(): confirmed.append(int(x[4]))
        if x[5].isdigit(): deaths.append(int(x[5]))
        if len(x)>6 and x[6].isdigit(): recovered.append(int(x[6]))

  smooth(confirmed,smoothlevel)
  smooth(deaths,smoothlevel)
  smooth(recovered,smoothlevel)
  for (c,d,r) in zip(confirmed,deaths,recovered):
    active.append(c-d-r)
  prev=0
  for c in confirmed:
    newc.append(c-prev);prev=c
  prev=0
  for d in deaths:
    newd.append(d-prev);prev=d

  for (i,x) in enumerate(confirmed):
    if x>=thr: break
  else: i=len(confirmed)

  # If smooth then remove last entry because it should be augmented by an unknown pullback from the future
  n=len(dates)-smoothlevel

  return (dates[i:n], confirmed[i:n], deaths[i:n], recovered[i:n], active[i:n], newc[i:n], newd[i:n])

def getallsimpledata(source="worldometer"):
  """
  Returns map (dict) from country to list of (date, confirmed cases, deaths, recovered, active), for all available countries
  If returned numerical values are of int type then valid, else will be '?'.
  """
  
  equivnames={}
  with open("countrynames") as fp:
    r=csv.reader(fp)
    for x in r:
      if len(x)==0 or x[0][:1]=='#': continue
      for y in x[1:]: equivnames[y]=x[0]

  results={}
  with open(source+".csv") as f:
    r=csv.reader(f)
    first=1
    dates=[]
    confirmed=[]
    deaths=[]
    recovered=[]
    active=[]
    for x in r:
      if first: first=0;continue
      country=equivnames.get(x[1],x[1])
      if country not in results: results[country]=[]
      date=x[0]
      confirmed=deaths=recovered=active='?'
      if x[4].isdigit(): confirmed=int(x[4])
      if x[5].isdigit(): deaths=int(x[5])
      if len(x)>6 and x[6].isdigit(): recovered=int(x[6])
      if confirmed!='?' and recovered!='?': active=confirmed-recovered
      results[country].append((date,confirmed,deaths,recovered,active))
  
  return results

# Adjust vv[] to be as smooth as possible by shifting readings back at most 1 day.  I.e.,
# we assume that there was a delay of 1 day in reporting some values and try to
# make the increment vector log(vv[i]/vv[i-1]) as near concave as possible, fixing vv[i0] (first non-zero entry) and vv[n-1].
def smooth(vv,lev):
  # Enforce monotonicity (covers bug in Norway data for 2020-03-07)
  prev=-1e30
  for i in range(len(vv)):
    if vv[i]<prev: vv[i]=prev
    else: prev=vv[i]
  if lev==0: return
  deb=False
  n=len(vv)
  for i0 in range(n):
    if vv[i0]>0: break
  else: return
  if i0==n-1: return
  dv=[0]*i0
  prev=1
  for i in range(i0,n): dv.append(log(vv[i]/prev));prev=vv[i]
  sh=dv[:]# To remember amount remaining available to be shifted back by 1 day (so we don't chain-shift more than 1 day)
  if deb: print("i0 =",i0,", vv =",vv)
  mxd=-1e9
  for i in range(i0+1,n):
    if dv[i]>mxd: mxd=dv[i];imx=i
  if deb: print("imx =",imx)
  chg=1e9;it=0
  while chg>1e-3 and it<1000:
    if deb:
      with open("tempd/out%d"%it,"w") as fp:
        for i in range(i0): print("%12g %12g"%(0,0),file=fp)
        t=0
        for i in range(i0,n): t+=dv[i];print("%12g %12g"%(dv[i],exp(t)),file=fp)
      for d in dv: print(" %12g"%d,end="")
      print(" %g"%chg)
    chg=0
    # Increasing bit
    for i in range(max(i0,1),imx-1):
      # Modify i,i+1 to improve i-1,i,i+1 (possibly worsening i,i+1,i+2 but that tends to have larger values to improve itself with later)
      s=max(min((dv[i-1]+dv[i+1]-2*dv[i])/3,(dv[i+1]-dv[i])/2,sh[i+1]),0)
      dv[i]+=s;dv[i+1]-=s;sh[i+1]-=s
      chg+=s
    # Decreasing bit
    for i in range(n-2,imx,-1):
      # Modify i,i+1 to improve i-1,i,i+1 (possibly worsening i-2,i-1,i but that tends to have larger values to improve itself with later)
      s=max(min((dv[i-1]+dv[i+1]-2*dv[i])/3,(dv[i-1]-dv[i])/2,sh[i+1]),0)
      dv[i]+=s;dv[i+1]-=s;sh[i+1]-=s
      chg+=s
    it+=1
    if deb: print(it,chg)
  t=0
  if deb: print("vv")
  for i in range(i0,n):
    t+=dv[i];vv[i]=exp(t)
    if deb: print(" %12g"%vv[i])
  if deb: print()


def getcountrylist(source="worldometer"):
  """
  Gets list of all countries for which there is information
  """
  
  equivnames={}
  with open("countrynames") as fp:
    r=csv.reader(fp)
    for x in r:
      if len(x)==0 or x[0][:1]=='#': continue
      for y in x[1:]: equivnames[y]=x[0]

  countries=set()
  with open(source+".csv") as f:
    r=csv.reader(f)
    first=1
    for x in r:
      if first: first=0;continue
      countries.add(equivnames.get(x[1],x[1]))

  return sorted(list(countries))

