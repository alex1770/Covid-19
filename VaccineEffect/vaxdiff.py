import time,calendar,os,json,sys,datetime,csv,random
from math import log,sqrt
from requests import get
from subprocess import Popen,PIPE

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def getapidata(req):
  url='https://api.coronavirus.data.gov.uk/v1/data?'
  response = get(url+req, timeout=10)
  if not response.ok:
    raise RuntimeError(f'Request failed: { response.text }')
  if response.status_code!=200: return None
  date=time.strftime('%Y-%m-%d',time.strptime(response.headers['Last-Modified'],'%a, %d %b %Y %H:%M:%S %Z'))# Not currently used
  return response.json()['data']

startdate='2020-11-01'

ltladatadir='ltladatadir'
stpdatadir='stpdatadir'
startday=datetoday(startdate)
nowdate=datetime.datetime.utcnow().strftime('%Y-%m-%d')
endday=datetoday(nowdate)+1

# Fetch and save new case data
os.makedirs(ltladatadir,exist_ok=True)
for day in range(startday,endday):
  date=daytodate(day)
  fn=os.path.join(ltladatadir,date)
  if os.path.isfile(fn): continue
  req='filters=areaType=ltla;date='+date+'&structure={"date":"date","LADCD":"areaCode","LADNM":"areaName","cases":"newCasesBySpecimenDateAgeDemographics"}'
  casedata=getapidata(req)
  if casedata==None: endday=day;break
  # Convert to a map from {LTLA location code} --> list of (min age, max age, case count), where [min age, max age) partition [0,150) in half-open convention
  data={}
  for x in casedata:
    data[x['LADCD']]=[]
    d=data[x['LADCD']]
    # Keep age ranges {xx_yy where xx is a multiple of 5 and <90, and yy=xx+4} and '90+'; discard others
    for y in x['cases']:
      a=y['age']
      if a=='90+':
        amin=90;amax=150
      elif '_' in a:
        amin=int(a[:2]);amax=int(a[-2:])
        if amax!=amin+4: continue
        amax+=1
      else: continue
      d.append((amin,amax,y['cases']))
  with open(fn,'w') as fp:
    json.dump(data,fp)
    print("Wrote",fn)

caseagebands=[(amin,amin+5 if amin<90 else 150) for amin in range(0,95,5)]
enddate=daytodate(endday)
print("Loading case counts from",startdate,"-",enddate,"(excluding end date)")

# Read LTLA->STP mapping
LTLAtoSTP={}
with open('LTLAtoSTP20.csv','r') as fp:
  r=csv.reader(fp)
  headings=next(r)
  (i0,i1)=(headings.index('LAD20CD'),headings.index('STP20NM'))
  for x in r:
    LTLAtoSTP[x[i0]]=x[i1]

# Need to modify conversion to use LAD19CD, because government api uses LAD19 but vaccination database (NIMS) uses STP20
for x in ['E07000004', 'E07000005', 'E07000006', 'E07000007']:
  LTLAtoSTP[x]=LTLAtoSTP['E06000060']

# Convert new LTLA case data to STP and save    
os.makedirs(stpdatadir,exist_ok=True)
for day in range(startday,endday):
  date=daytodate(day)
  fnstp=os.path.join(stpdatadir,date)
  if os.path.isfile(fnstp): continue
  fnltla=os.path.join(ltladatadir,date)
  with open(fnltla,'r') as fp:
    ltladata=json.load(fp)
  stpdata={}
  for ltla in ltladata:
    cases=ltladata[ltla]
    stp=LTLAtoSTP[ltla]
    if stp not in stpdata: stpdata[stp]=[]
    for (i,c) in enumerate(cases):
      if len(stpdata[stp])==i: stpdata[stp].append(cases[i])
      else: stpdata[stp][i][2]+=cases[i][2]
  with open(fnstp,'w') as fp:
    json.dump(stpdata,fp)
    print("Wrote",fnstp)

# Load all STP case data
# stpdat[day number from startdate][STP name][age band number] = new case count for that day, place, age band
stpcases=[]
for day in range(startday,endday):
  date=daytodate(day)
  fnstp=os.path.join(stpdatadir,date)
  fnstp=os.path.join(stpdatadir,date)
  with open(fnstp,'r') as fp:
    stpcases.append(json.load(fp))

# Load STP population size data
stppop={}
with open('STPpopulations.csv','r') as fp:
  r=csv.reader(fp)
  headings=next(r)
  for x in r:
    stppop[x[1]]=[]
    for (h,v) in zip(headings[2:],x[2:]):
      if '-' in h:
        z=h.split('-')
        amin,amax=int(z[0]),int(z[1])+1
      elif '+' in h:
        amin,amax=int(h[:-1]),150
      else: assert 0
      stppop[x[1]].append((amin,amax,int(v)))

# Check STP names are compatible
assert set(stppop.keys())==set(LTLAtoSTP.values())
STPNM=sorted(list(set(LTLAtoSTP.values())))

# Read STP -> NHS region mapping
STPtoNHS={}
with open('STP20toNHSER20.csv','r') as fp:
  r=csv.reader(fp)
  headings=next(r)
  (i0,i1)=(headings.index('STP20NM'),headings.index('NHSER20NM'))
  for x in r:
    STPtoNHS[x[i0]]=x[i1]
NHSER=sorted(list(set(STPtoNHS.values())))

# Load vaccination data
# Weekly data from https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/
# Note that publication date is 4 days after date of last data. Filename uses the latter date.
# Hard code this particular week until we decide we want to do something else.
# vaxnum[STP name] = list of (agemin, agemax, count)
vaxnum={}
vaxagebands=[(0,80),(80,150)]
vaxdate='2021-01-17'
with open(vaxdate+'-vaxdata.csv','r') as fp:
  r=csv.reader(fp)
  headings=next(r)
  # For now assume headings are: Region, STP name, 1st dose 0-80, 1st dose 80+, 2nd dose 0-80, 2nd dose 80+
  for x in r:
    vaxnum[x[1]]=[(0,80,int(x[2])),(80,150,int(x[3]))]
print("Using vaccine data from",vaxdate)

# Check STP names are compatible
assert set(vaxnum)==set(stppop.keys())

# Smooth cases with lagged 7 day average
newstpcases=[]
newstartday=startday+6
for day in range(newstartday,endday):
  newstpcases.append({})
  for stp in STPNM:
    newstpcases[-1][stp]=[(amin,amax,sum(stpcases[d-startday][stp][i][2] for d in range(day-6,day+1))/7) for (i,(amin,amax)) in enumerate(caseagebands)]
stpcases=newstpcases
startday=newstartday
startdate=daytodate(startday)

from scipy.optimize import minimize
import numpy as np

# Add up values in l corresponding to intervals contained in [rmin,rmax)
def totrange(l,rmin,rmax):
  return sum(x[2] for x in l if x[0]>=rmin and x[1]<=rmax)

def sc(xx):
  eff,ff=xx[0],xx[1:]
  err=0
  for stp in STPNM:
    for (amin,amax),f in zip(vaxagebands,ff):
      #if amin==0: continue
      t0=totrange(stpcases[day0-startday][stp],amin,amax)
      t1=totrange(stpcases[day1-startday][stp],amin,amax)
      vax=totrange(vaxnum[stp],amin,amax)
      pop=totrange(stppop[stp],amin,amax)
      v=vax/pop
      #v=randvax[(stp,amin)]
      rawpred=t0
      actual=t1
      pred=f*(1-v*eff)*rawpred
      #print("%-20s   %3d %3d       %5.3f     %6d  %6d  %6.1f"%(stp[:20].replace(' ','_'),amin,amax,v,t0,actual,pred))
      err+=(log(pred/actual))**2
      #err-=-pred+actual*log(pred)
  return err

randvax={}
#for stp in STPNM:
#  for (amin,amax) in vaxagebands:
#    randvax[(stp,amin)]=0.3+0.5*random.random()
for nhs in NHSER:
  for (amin,amax) in vaxagebands:
    r=0.3+0.5*random.random()
    for stp in STPNM:
      if STPtoNHS[stp]==nhs: randvax[(stp,amin)]=r

def sc2(xx):
  eff,f=xx[0],xx[1]
  err=0
  for stp in STPNM:
    casethen=[]
    casenow=[]
    vaxrate=[]
    for (amin,amax) in vaxagebands:
      t0=totrange(stpcases[day0-startday][stp],amin,amax)
      t1=totrange(stpcases[day1-startday][stp],amin,amax)
      vax=totrange(vaxnum[stp],amin,amax)
      pop=totrange(stppop[stp],amin,amax)
      casethen.append(t0)
      casenow.append(t1)
      vaxrate.append(vax/pop)
      #vaxrate.append(randvax[(stp,amin)])
    v=vaxrate[1]
    # Cross-ratio + Poisson model
    pred=casenow[0]/casethen[0]*casethen[1]*f*(1-v*eff)
    actual=casenow[1]
    err-=-pred+actual*log(pred)# I don't think Poisson model is very good here
    #err+=log(pred/actual)**2
  return err

if 0:
  day1=endday-1
  print("Measuring changes in case counts up to",daytodate(day1),"(7 day lagged average)")
  print()
  best=(-1,)
  print("From date    Efficacy coefficient")
  for day0 in range(datetoday(vaxdate)-10,datetoday(vaxdate)+15):
    res=minimize(sc,[0.5, 0.5, 0.5],method="SLSQP",bounds=[(0,0.99),(1e-6,1),(1e-6,1)],options={"maxiter":1000})
    #res=minimize(sc2,[0.5, 1],method="SLSQP",bounds=[(0,0.99),(1e-6,10),],options={"maxiter":1000})
    if not res.success: raise RuntimeError(res.message)
    xx=res.x
    eff,ff=xx[0],xx[1:]
    if eff>best[0]: best=(eff,day0,ff)
    print(daytodate(day0),"  %5.3f"%eff)
  print()

if 0:
  day1s=list(reversed(range(endday-1,datetoday('2020-11-25'),-4)))
  print("From date  "+' '.join(daytodate(x) for x in day1s))
  for day0 in reversed(range(endday-1,datetoday(vaxdate)-50,-2)):
    print(daytodate(day0),end='')
    for day1 in day1s:
      if day1<=day0: print("           ",end='');continue
      res=minimize(sc2,[0.5, 1],method="SLSQP",bounds=[(-9,0.99),(1e-6,10),],options={"maxiter":1000})
      if not res.success: raise RuntimeError(res.message)
      l0=-res.fun
      eff,f=res.x
      if 1:
        eff1=eff
      else:
        beta=0.8
        res=minimize(sc2,[eff*beta, f],method="SLSQP",bounds=[(eff*beta,eff*beta),(1e-6,10)],options={"maxiter":1000}) 
        if not res.success: raise RuntimeError(res.message)
        l1=-res.fun
        if l0<=l1: eff1=0
        else:
          eff1=eff-((1-beta)*eff)/sqrt((l0-l1)/20)
          if eff*eff1<0: eff1=0
      print("     %6.3f"%(eff1),end='')
    print()
  print()

if 1:
  day1s=list(reversed(range(endday-1,datetoday('2020-11-25'),-4)))
  print("From date  "+' '.join(daytodate(x) for x in day1s))
  for day0 in reversed(range(endday-1,datetoday(vaxdate)-50,-2)):
    print(daytodate(day0),end='')
    for day1 in day1s:
      if day1<=day0: print("           ",end='');continue
      res=minimize(sc,[0.5, 0.5, 0.5],method="SLSQP",bounds=[(0,0.99),(1e-6,10),(1e-6,10)],options={"maxiter":1000})
      if not res.success: raise RuntimeError(res.message)
      l0=-res.fun
      eff,ff=res.x[0],res.x[1:]
      if 0:
        eff1=eff
      else:
        beta=0.8
        res=minimize(sc,[eff*beta, ff[0], ff[1]],method="SLSQP",bounds=[(eff*beta,eff*beta),(1e-6,10),(1e-6,10)],options={"maxiter":1000}) 
        if not res.success: raise RuntimeError(res.message)
        l1=-res.fun
        if l0<=l1: eff1=0
        else:
          eff1=eff-((1-beta)*eff)/sqrt((l0-l1)/2)
          if eff*eff1<0: eff1=0
      print("     %6.3f"%(eff1),end='')
    print()
  print()

if 0:
  (eff,day0,ff)=best
  print("Strongest signal in changes in case counts is from",daytodate(day0),"-",daytodate(day1),"(in both cases smoothed with 7 day lagging average)")
  print("Coefficient of efficacy : %5.3f"%eff)
  for (amin,amax),f in zip(vaxagebands,ff):
    print("Case count factor for age band %3d - %3d : %5.3f"%(amin,amax,f))
  print()

if 0:
  #day0=datetoday('2020-12-12')
  #day1=datetoday('2021-01-06')
  res=minimize(sc2,[0.5, 1],method="SLSQP",bounds=[(-9,0.99),(1e-6,10),],options={"maxiter":1000})
  if not res.success: raise RuntimeError(res.message)
  xx=res.x
  eff,ff=xx[0],xx[1:]
  for stp in STPNM:
    for (amin,amax),f in zip(vaxagebands,ff):
      t0=totrange(stpcases[day0-startday][stp],amin,amax)
      t1=totrange(stpcases[day1-startday][stp],amin,amax)
      vax=totrange(vaxnum[stp],amin,amax)
      pop=totrange(stppop[stp],amin,amax)
      v=vax/pop
      print("%5d  %6.1f  %5.3f  %5d  %6.3f  %d"%(amin,t0*f,v,t1,log(t1/(t0*f)),NHSER.index(STPtoNHS[stp])))

if 0:
  day0=datetoday('2021-01-24')
  day1=endday-1
  #day0=datetoday('2020-12-12')
  #day1=datetoday('2021-01-06')
  res=minimize(sc2,[0.5, 1],method="SLSQP",bounds=[(-9,0.99),(1e-6,10),],options={"maxiter":1000})
  #res=minimize(sc2,[0, 1],method="SLSQP",bounds=[(0,0),(1e-6,10),],options={"maxiter":1000})
  if not res.success: raise RuntimeError(res.message)
  xx=res.x
  eff,f=xx
  print("eff =",eff)
  print("f =",f)
  totl=0
  fp=open('temp','w')
  for stp in STPNM:
    casethen=[]
    casenow=[]
    vaxrate=[]
    print("%25s"%(STPtoNHS[stp]),end='')
    for (amin,amax) in vaxagebands:
      t0=totrange(stpcases[day0-startday][stp],amin,amax)
      t1=totrange(stpcases[day1-startday][stp],amin,amax)
      vax=totrange(vaxnum[stp],amin,amax)
      pop=totrange(stppop[stp],amin,amax)
      casethen.append(t0)
      casenow.append(t1)
      vaxrate.append(vax/pop)
      print("   %6d  %5.3f"%(vax,vax/pop),end='')
    print()
    v=vaxrate[1]
    rawpred=casenow[0]/casethen[0]*casethen[1]
    pred=rawpred*f*(1-v*eff)
    actual=casenow[1]
    #print("%6.1f  %6.1f   %d"%(pred,actual,NHSER.index(STPtoNHS[stp])))
    l=-pred+actual*log(pred)
    totl+=l
    print("%6.3f  %6.3f   %6.3f   %12g %12g %12g %d   %12g"%(v,actual/rawpred,actual/pred,actual,rawpred,pred,NHSER.index(STPtoNHS[stp]),l),file=fp)
  print("Total likelihood",totl)
  fp.close()

sys.exit(0)

with open('stpcrossratios','w') as fp:
  for stp in STPNM:
    casethen=[]
    casenow=[]
    vaxrate=[]
    for (amin,amax),f in zip(vaxagebands,ff):
      t0=totrange(stpcases[day0-startday][stp],amin,amax)
      t1=totrange(stpcases[day1-startday][stp],amin,amax)
      vax=totrange(vaxnum[stp],amin,amax)
      pop=totrange(stppop[stp],amin,amax)
      casethen.append(t0)
      casenow.append(t1)
      vaxrate.append(vax/pop)
    print('%5.3f   %6.3f  "%s"'%(vaxrate[1],log((casenow[1]/casethen[1])/(casenow[0]/casethen[0])),stp),file=fp)

with open('nhscrossratios','w') as fp:
  for reg in NHSER:
    casethen=[]
    casenow=[]
    vaxrate=[]
    for (amin,amax),f in zip(vaxagebands,ff):
      t0=t1=vax=pop=0
      for stp in STPNM:
        if STPtoNHS[stp]!=reg: continue
        t0+=totrange(stpcases[day0-startday][stp],amin,amax)
        t1+=totrange(stpcases[day1-startday][stp],amin,amax)
        vax+=totrange(vaxnum[stp],amin,amax)
        pop+=totrange(stppop[stp],amin,amax)
      casethen.append(t0)
      casenow.append(t1)
      vaxrate.append(vax/pop)
    print('%5.3f   %6.3f  "%s"'%(vaxrate[1],log((casenow[1]/casethen[1])/(casenow[0]/casethen[0])),reg),file=fp)

# gnuplot> plot "stpcrossratios" u 1:2:(strcol(3)[1:12]) with labels point pt 7 offset char 0,-1.0
# gnuplot> plot "nhscrossratios" u 1:2:3 with labels point pt 7 offset char 0,-1.2
