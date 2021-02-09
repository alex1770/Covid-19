import time,calendar,os,json,sys,datetime,csv
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

startdate='2020-12-01'

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

enddate=daytodate(endday)
print("Using dates",startdate,"-",enddate,"(excluding end date)")

# Read LTLA->STP mapping
LTLAtoSTP={}
with open('LTLAtoSTP20','r') as fp:
  r=csv.reader(fp)
  headings=next(r)
  (i0,i1)=(headings.index('LAD20CD'),headings.index('STP20NM'))
  for x in r:
    LTLAtoSTP[x[i0]]=x[i1]

# Need to modify conversion to use LAD19CD, because vaccination database (NIMS) uses a mixture of LAD19 and STP20
for x in ['E07000004', 'E07000005', 'E07000006', 'E07000007']:
  LTLAtoSTP[x]=LTLAtoSTP['E06000060']

# Convert LTLA case data to STP and save    
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

# Load STP population data
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

# Load vaccination data
# Weekly data from https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/
# Note that publication date is 4 days after date of last data. Filename uses the latter date.
# Hard code this particular week for the moment.
# vaxnum[STP name] = list of (agemin, agemax, count)
vaxnum={}
with open('2021-01-17-vaxdata.csv','r') as fp:
  r=csv.reader(fp)
  headings=next(r)
  # For now assume headings are: Region, STP name, 1st dose 0-80, 1st dose 80+, 2nd dose 0-80, 2nd dose 80+
  for x in r:
    vaxnum[x[1]]=[(0,80,int(x[2])),(80,150,int(x[3]))]

