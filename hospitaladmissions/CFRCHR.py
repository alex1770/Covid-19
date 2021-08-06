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
  x=x.replace('_to_','_').replace(' to ','_')# cater for 65_to_69, 65_69, "65 to 69" formats
  aa=[int(y) for y in x.split("_")]
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

# Interpolate cases to hospages bands, (0,6) (6,18) (18,65) (65,85) (85,150), based on school assumption: I.e., that 0-5, 5-18, 18+ groups are similar
print("Cases -> Hospitalisations using dashboard hosp figures")
ages=hospages+[(0,150)]
for h in hosps:
  h[(0,150)]=sum(h[x] for x in hospages)
for c in cases:
  c[(0,6)]=c[(0,5)]+c[(5,10)]/5
  c[(6,18)]=c[(5,10)]*4/5+c[(10,15)]*8/5
  c[(18,65)]=c[(20,25)]*2/5+c[(20,25)]+c[(25,30)]+c[(30,35)]+c[(35,40)]+c[(40,45)]+c[(45,50)]+c[(50,55)]+c[(55,60)]+c[(60,65)]
  for a in [(65,85),(85,150)]:
    c[a]=sum(c[x] for x in caseages if x[0]>=a[0] and x[1]<=a[1])
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

#################################################################################################################################
# Cases -> Hospitalisations using hosp figures from National flu and COVID-19 surveillance report (Fig.40 SARIWatch-hospagegrp) #
# https://www.gov.uk/government/statistics/national-flu-and-covid-19-surveillance-reports                                       #
#################################################################################################################################

print()
print("Cases -> Hospitalisations using National flu and COVID-19 surveillance report figures")
h2=loadcsv("hosps_byage.csv")
conv={}
for k in h2.keys():
  try:
    a=parseage(k)
    conv[k]=a
  except:
    pass

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

hospages2=sorted(list(conv.values()))
ages=hospages2+[(0,150)]
popages=list(pop)
for c in cases:
  for a in ages:
    c[a]=sum(c[x] for x in caseages if x[0]>=a[0] and x[1]<=a[1])
for p in [pop]:
  for a in ages:
    p[a]=sum(p[x] for x in popages   if x[0]>=a[0] and x[1]<=a[1])

n=len(h2['w/e'])
hosps2=[]
for i in range(n):
  h={}
  h['date']=h2['w/e'][i]
  for k in conv:
    h[conv[k]]=int(h2[k][i]*p[conv[k]]/1e5+.5)
  hosps2.append(h)
    
for h in hosps2:
  for a in ages:
    h[a]=sum(h[x] for x in hospages2 if x[0]>=a[0] and x[1]<=a[1])

def inc(date): return daytodate(datetoday(date)+1)
    
# Dates of hosps2[-i] are those of cases[-i*7+offset+{0,...,6}] (and corresponds to precases[-i])
offset=datetoday(hosps2[-1]['date'])-datetoday(cases[-1]['date'])
precases=[{} for i in range(len(hosps2))]
ave=7 # perforce, since this hosp data is grouped in weeks
minlag=4
maxlag=8
back=180
for n in range(-back//7-2,0):
  for a in ages:
    t=0
    for i in range(-maxlag,-minlag):
      for j in range(7):
        t+=cases[n*7+offset+j+i][a]
    precases[n][a]=t/(maxlag-minlag)

n0=-6
for n in range(-back//7,0):
  print(inc(hosps2[n-1]['date']),'-',hosps2[n]['date'],end='')
  for a in ages:
    h=hosps2[n][a]
    print(" %6d"%h,end='')
  print()
  print(inc(hosps2[n-1]['date']),'-',hosps2[n]['date'],end='')
  for a in ages:
    p=precases[n][a]
    print(" %6d"%p,end='')
  print()
  print(inc(hosps2[n-1]['date']),'-',hosps2[n]['date'],end='')
  for a in ages:
    h=hosps2[n][a]
    p=precases[n][a]
    print(" %6.2f"%(h/p*100),end='')
  print()
  if n>=n0:
    print(inc(hosps2[n-1]['date']),'-',hosps2[n]['date'],end='')
    for a in ages:
      h=hosps2[n][a]
      p=precases[n][a]
      print(" %6.2f"%(h/p/(hosps2[n0][a]/precases[n0][a])),end='')
    print()
  print()
print("                       ",end='')
for a in ages: print(" %6s"%("%d-%d"%a),end='')
print()

