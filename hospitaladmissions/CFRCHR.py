import time,calendar,os,json,sys,datetime
from requests import get
from math import sqrt,log,exp
from stuff import *

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v1/data?'
  response = get(url+req, timeout=10)
  if not response.ok:
    raise RuntimeError(f'Request failed: { response.text }')
  date=time.strftime('%Y-%m-%d',time.strptime(response.headers['Last-Modified'],'%a, %d %b %Y %H:%M:%S %Z'))# Not currently used
  data=response.json()['data']
  
  # Convert from list form to dictionary keyed by age
  day=datetoday(data[0]['date'])
  n=1
  while n<len(data) and datetoday(data[n]['date'])==day-n: n+=1# Find maximal contiguous date range
  data1=[]
  for i in range(n-1,-1,-1):
    d=data[i]
    e={'date':d['date']}
    for x in d:
      if x!='date':
        for y in d[x]:
          if 'value' in y: val=y['value']
          else: val=y['deaths']
          e[y['age']]=e.get(y['age'],0)+val
    data1.append(e)
  return data1

req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"newDeaths28DaysByDeathDateAgeDemographics"}'; mortdata=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"cumAdmissionsByAge"}';                        hospdata=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","male":"maleCases"}';                                 malecases=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","female":"femaleCases"}';                             femalecases=get_data(req)

# Interpolate to get v[(a,b)] for (a,b) in l. Assumed that p is a partition = [(a0,a1), (a1,a2), (a2,a3), ...]
def interpolate(v,p,l):
  for (a,b) in l:
    t=0
    for (c,d) in p:
      if b<=c or d<=a: continue
      t+=(min(b,d)-max(a,c))/(d-c)*v[(c,d)]
    v[(a,b)]=t
  
casedata=[]
for (m,f) in zip(malecases,femalecases):
  d={'date': m['date']}
  assert m['date']==f['date']
  for s in [m,f]:
    for x in s:
      if x!='date': d[x]=d.get(x,0)+s[x]
  casedata.append(d)

def getdiff(data):
  n=len(data)
  newdata=[]
  for i in range(1,n):
    l={'date':data[i]['date']}
    for age in data[i]:
      if age!='date': l[age]=data[i][age]-data[i-1].get(age,0)
    newdata.append(l)
  return newdata

newhosps=getdiff(hospdata)
newcases=getdiff(casedata)
newmcases=getdiff(malecases)
newfcases=getdiff(femalecases)
discard=2# Last entries are unreliable.
newcases=newcases[:-discard]
newmcases=newmcases[:-discard]
newfcases=newfcases[:-discard]
                  
# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','-').replace(' to ','-').replace('_','-')# cater for 65_to_69, 65_69, "65 to 69", 15-19 formats
  aa=[int(y) for y in x.split("-")]
  return (aa[0],aa[1]+1)

# Convert dictionary from using '15_19' (etc) format to (15,20) format
# At the same time remove age brackets such as '60+' and '00_59' that strictly contain other age brackets, so avoiding overcounting
# Return list of ages
def convertages(dd):
  ages0=[(x,parseage(x)) for x in dd[-1] if x!='date']
  ages1={}
  for (x,(a,b)) in ages0:
    for (y,(c,d)) in ages0:
      if c>=a and d<=b and (c>a or d<b): break
    else: ages1[x]=(a,b)
  ee=[]
  for d in dd:
    e={}
    e['date']=d['date']
    for x in ages1:
      e[ages1[x]]=d.get(x,0)
    ee.append(e)
  ages2=sorted(ages1.values())
  return (ee,ages2)

hosps,hospages=convertages(newhosps)
cases,caseages=convertages(newcases)
deaths,deathages=convertages(mortdata)
fcases,_=convertages(newfcases)
mcases,_=convertages(newmcases)

assert caseages==deathages
ages=[(0,40)]+[(i,i+10) for i in range(40,90,10)]+[(90,150),(0,150)]
for c in cases:
  for a in ages:
    c[a]=sum(c[x] for x in caseages if x[0]>=a[0] and x[1]<=a[1])
for d in deaths:
  for a in ages:
    d[a]=sum(d[x] for x in deathages if x[0]>=a[0] and x[1]<=a[1])

###################
# Cases -> Deaths #
###################

print("Cases -> Deaths")
offset=datetoday(deaths[-1]['date'])-datetoday(cases[-1]['date'])
precases=[{} for i in range(len(deaths))]
ave=7
minlag=15
maxlag=21
back=180
for n in range(-back-ave,0):
  for a in ages:
    t=0
    for i in range(-maxlag,-minlag):
      t+=cases[n+i+offset][a]
    precases[n][a]=t/(maxlag-minlag)

for n in range(-back,0):
  print(deaths[n-(ave-1)]['date'],'-',deaths[n]['date'],end='')
  for a in ages:
    d=sum(deaths[n-i][a] for i in range(ave))
    p=sum(precases[n-i][a] for i in range(ave))
    print(" %6.2f"%(d/p*100),end='')
  print()
print("                       ",end='')
for a in ages: print(" %6s"%("%d-%d"%a),end='')
print()

print()
for n in range(-back,-1+back-1,back-1):
  print(deaths[n-(ave-1)]['date'],'-',deaths[n]['date'],end='')
  for a in ages:
    d=sum(deaths[n-i][a] for i in range(ave))
    print(" %6d"%d,end='')
  print()
  print(deaths[n-(ave-1)]['date'],'-',deaths[n]['date'],end='')
  for a in ages:
    p=sum(precases[n-i][a] for i in range(ave))
    print(" %6.0f"%p,end='')
  print()
print("                       ",end='')
for a in ages: print(" %6s"%("%d-%d"%a),end='')
print()
print()

##########################################################
# Cases -> Hospitalisations using dashboard hosp figures #
##########################################################

print("Cases -> Hospitalisations using dashboard hosp figures")
ages=hospages+[(0,150)]
for h in hosps:
  h[(0,150)]=sum(h[x] for x in hospages)

for c in cases:
  interpolate(c,caseages,[(0,6),(20,65),(25,35),(35,45),(45,55),(55,65),(65,75),(75,85),(65,85),(85,150)])
  # Special hand-picked interpolation for these bands because want to assume that school children (5-18) are similar to each other
  # and student age (18-25) are also similar to each other:
  c[(6,18)]=c[(5,10)]*4/5+c[(10,15)]+c[(10,15)]*3/5
  c[(18,25)]=c[(20,25)]*2/5+c[(20,25)]
  c[(18,65)]=c[(20,25)]*2/5+c[(20,65)]
  #c[(0,150)]=sum(c[x] for x in caseages)

offset=datetoday(hosps[-1]['date'])-datetoday(cases[-1]['date'])
precases=[{} for i in range(len(hosps))]
ave=7
minlag=4
maxlag=8
back=180
for n in range(-back-ave,0):
  for a in ages:
    t=0
    for i in range(-maxlag,-minlag):
      t+=cases[n+i+offset][a]
    precases[n][a]=t/(maxlag-minlag)

for n in range(-back,0):
  if 0:
    print()
    print(hosps[n-(ave-1)]['date'],'-',hosps[n]['date'],end='')
    for a in ages:
      h=sum(hosps[n-i][a] for i in range(ave))
      print(" %6d"%h,end='')
    print()
    print(hosps[n-(ave-1)]['date'],'-',hosps[n]['date'],end='')
    for a in ages:
      p=sum(precases[n-i][a] for i in range(ave))
      print(" %6d"%p,end='')
    print()
  print(hosps[n-(ave-1)]['date'],'-',hosps[n]['date'],end='')
  for a in ages:
    h=sum(hosps[n-i][a] for i in range(ave))
    p=sum(precases[n-i][a] for i in range(ave))
    print(" %6.2f"%(h/p*100),end='')
  print()
print("                       ",end='')
for a in ages: print(" %6s"%("%d-%d"%a),end='')
print()

#########################################################################################################################################
# Cases -> Hospitalisations using hosp figures from (changes each week/month?)                                                          #
# https://www.england.nhs.uk/statistics/wp-content/uploads/sites/2/2021/08/Covid-Publication-05-08-2021-Supplementary-Dataagebands.xlsx #
# contained in https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-hospital-activity/                                 #
#########################################################################################################################################

print()
print("Cases -> Hospitalisations, using NHS admissions analysis by age group")
h2=loadcsv("Englandadmissionsbyfineage.csv")
age2ind={}
for (i,s) in enumerate(h2['Measure']):
  t=s.split()[-1]
  if t[0].isdigit():
    age2ind[parseage(t)]=i
hospages2=sorted(list(age2ind))
ages=hospages2+[(0,150)]

hdays=sorted(datetoday(x) for x in h2.keys() if x!='Measure')
nhdays=len(hdays)
minhday=hdays[0]
assert hdays==list(range(minhday,minhday+nhdays))# Checking contiguous
hospcases2=[{} for i in range(nhdays)]
for x in h2:
  if x!='Measure':
    i=datetoday(x)-minhday
    for a in hospages2:
      hospcases2[i][a]=h2[x][age2ind[a]]
    interpolate(hospcases2[i],hospages2,[(0,150)])

pop={     (0,5):   3937768,
         (5,10):   4116901,
        (10,15):   3920483,
        (15,20):   3685333,
        (20,25):   4097454,
        (25,30):   4497206,
        (30,35):   4695232,
        (35,40):   4554064,
        (40,45):   4289457,
        (45,50):   4325007,
        (50,55):   4654729,
        (55,60):   4501186,
        (60,65):   3852250,
        (65,70):   3394375,
        (70,75):   3358992,
        (75,80):   2393054,
        (80,85):   1717798,
        (85,90):   1070495,
        (90,150):   646294}

popages=sorted(list(pop))
for p in [pop]:
  for a in ages:
    p[a]=sum(p[x] for x in popages   if x[0]>=a[0] and x[1]<=a[1])
interpolate(pop,popages,ages)

# Dates of hosps2[-i] are those of cases[-i*7+offset+{0,...,6}] (and corresponds to precases[-i])
offset=minhday+nhdays-1-datetoday(cases[-1]['date'])
ave=1
minlag=4
maxlag=8
back=180
precases=[{} for i in range(back//ave+3)]
for n in range(-back//ave-2,0):
  for a in ages:
    t=0
    for i in range(-maxlag,-minlag):
      for j in range(ave):
        t+=cases[n*ave+offset+j+i][a]
    precases[n][a]=t/(maxlag-minlag)/ave

hospave=[{} for i in range(back//ave+1)]
for n in range(-back//ave,0):
  for a in ages:
    t=0
    for j in range(ave):
      t+=hospcases2[n*ave+j][a]
    hospave[n][a]=t/ave
    
for n in range(-back//ave,0):
  dates=daytodate(minhday+nhdays+n*ave)+' - '+daytodate(minhday+nhdays+(n+1)*ave-1)
  print(dates,end='')
  for a in ages:
    h=hospave[n][a]
    print(" %6d"%h,end='')
  print("   Hospitalisations")
  print(dates,end='')
  for a in ages:
    p=precases[n][a]
    print(" %6d"%p,end='')
  print("   Cases")
  print(dates,end='')
  for a in ages:
    h=hospave[n][a]
    p=precases[n][a]
    print(" %6.2f"%(h/p*100),end='')
  print("   Percent hospitalised")
  print()
print("                       ",end='')
for a in ages: print(" %6s"%("%d-%d"%a),end='')
print()

