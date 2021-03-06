import apitools,json,os,time,calendar
import numpy as np

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def getcasedata():
  cdir="casedata"
  os.makedirs(cdir,exist_ok=True)
  laststored=max((x for x in os.listdir(cdir) if x[:2]=="20"),default="0000-00-00")
  query_structure = {
    "date": "date",
    "name": "areaName",
    "code": "areaCode",
    "cases": "newCasesByPublishDate",
  }
  # Small test query to get the latest date
  data = apitools.get_paginated_dataset(["areaType=ltla","areaName=Barnet"],query_structure)
  lastavail=data[0]['date']
  if lastavail>laststored:
    data = apitools.get_paginated_dataset(["areaType=ltla"], query_structure)
    with open(os.path.join(cdir,lastavail),'w') as fp:
      json.dump(data,fp,indent=2)
  else:
    with open(os.path.join(cdir,laststored),'r') as fp:
      data=json.load(fp)
  return data

data=getcasedata()

with open('zoemapdata/2021-01-17','r') as fp: zm=json.load(fp)

zoepl=dict((pl['lad16cd'],pl['lad16nm']) for pl in zm.values())
offpl=dict((x['code'],x['name']) for x in data)

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
both=set((x['date'],x['code']) for x in data if x['date']>=mindate)
r=0
for date in sorted(list(dates)):
  for loc in list(locs):
    if (date,loc) not in both and loc in locs: locs.remove(loc);r+=1
if r>0: print("Removed %d locations that didn't have up to date case counts"%r)

locs=sorted(list(locs))
locind=dict((loc,i) for (i,loc) in enumerate(locs))

# Convert Zoe and Official data into convenient sequences
N=len(locs)
minday=datetoday(mindate)
n=datetoday(lastdate)-minday+1
zvals=np.zeros([2,N,n])# (predicted cases)/(respondents)*(population), (corrected covid positive)
cases=np.zeros([N,n])
zdates=np.zeros(n,dtype=bool)
for d in data:
  t=datetoday(d['date'])-minday
  if t>=0 and d['code'] in locind: cases[locind[d['code']],t]=d['cases']
for date in os.listdir('zoemapdata'):
  if date[:2]=='20' and date>=mindate and date<=lastdate:
    with open(os.path.join('zoemapdata',date),'r') as fp:
      zm=json.load(fp)
      t=datetoday(date)-minday
      zdates[t]=1
      for d in zm.values():
        if d['lad16cd'] in locind:
          i=locind[d['lad16cd']]
          zvals[0,i,t]=d['predicted_covid_positive_count']/d['respondent']*d['population']
          zvals[1,i,t]=d['corrected_covid_positive']
# Interpolate missing dates
for t in range(n):
  if not zdates[t]:
    for t0 in range(t-1,-1,-1):
      if zdates[t0]: break
    for t1 in range(t+1,n):
      if zdates[t1]: break
    zvals[:,:,t]=((t1-t)*zvals[:,:,t0]+(t-t0)*zvals[:,:,t1])/(t1-t0)

for delta in range(-40,41):
  if delta>0:  zvals1=zvals[:,:,:-delta]
  else:        zvals1=zvals[:,:,-delta:]
  if delta>=0: cases1=cases[:,delta:]
  else:        cases1=cases[:,:delta]
  n1=n-delta
  
  # Maybe scale first by a location-dependent, time-independent function
  for i in range(zvals1.shape[0]):
    lam=np.tensordot(zvals1[i,:,:],cases1)/np.tensordot(zvals1[i,:,:],zvals1[i,:,:])
  
    errmat=lam*zvals1[i,:,:]-cases1
    err=np.tensordot(errmat,errmat)/1e6/zvals1.shape[2]
    print("%3d  %d  %7.3f"%(delta,i,err))
  #print()
