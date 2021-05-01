# Estimate how many (non-Covid) deaths that would normally have occurred on a given day from the confirmed cases between day-27 and day
# (remembering of course that confirmed cases is less than actual cases by some factor)
import time,calendar,os,json,sys
from stuff import *
from requests import get
import numpy as np

def get_data0(req):
  url='https://api.coronavirus.data.gov.uk/v1/data?'
  response = get(url+req, timeout=10)
  if not response.ok:
    raise RuntimeError(f'Request failed: { response.text }')
  data=response.json()['data']
  day=datetoday(data[0]['date'])
  n=1
  while n<len(data) and datetoday(data[n]['date'])==day-n: n+=1# Find maximal contiguous date range
  return data[n-1::-1]

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v1/data?'
  response = get(url+req, timeout=10)
  if not response.ok:
    raise RuntimeError(f'Request failed: { response.text }')
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
req='filters=areaType=nation;areaName=england&structure={"date":"date","deaths":"newDeaths28DaysByDeathDate"}'; totdeaths=get_data0(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","male":"maleCases"}';              malecases=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","female":"femaleCases"}';          femalecases=get_data(req)

casedata=[]
for (m,f) in zip(malecases,femalecases):
  d={'date': m['date']}
  assert m['date']==f['date']
  for s in [m,f]:
    for x in s:
      if x!='date': d[x]=d.get(x,0)+s[x]
  casedata.append(d)

updatedate=casedata[-1]['date']
now=datetime.datetime.utcnow().strftime('%Y-%m-%d')

def getdiff(data):
  n=len(data)
  newdata=[]
  for i in range(1,n):
    l={'date':data[i]['date']}
    for age in data[i]:
      if age!='date': l[age]=data[i][age]-data[i-1].get(age,0)
    newdata.append(l)
  return newdata

newcases=getdiff(casedata)
newmcases=getdiff(malecases)
newfcases=getdiff(femalecases)
newcases=newcases[:-1]# Last entry seems particularly unreliable, I think because it using specimen date and there are biases with recent entries
newmcases=newmcases[:-1]
newfcases=newfcases[:-1]
                  
# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

# Convert (eg) (15,20) to "15 - 19"
def unparse(r):
  (a,b)=r
  if b==150: return "%d+"%a
  return "%d - %d"%(a,b)

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

mindate='2020-06-01'
cases,caseages=convertages(newcases)
deaths,deathages=convertages(mortdata)
fcases,_=convertages(newfcases)
mcases,_=convertages(newmcases)

ageses=[sorted([x for x in l[-1].keys() if x!='date']) for l in [deaths,mcases,fcases]]
assert ageses[0]==ageses[1]==ageses[2]
ages=ageses[0]

# Death rates by age and sex from https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables
deathrates=loadcsv("mortalitybyageandsex.csv")
adr={}# Average death rates by age band; weight the ages within the band according to likely survival
for (desc,dr) in [('m',deathrates['MaleDeathRate']),('f',deathrates['FemaleDeathRate'])]:
  adr[desc]={}
  for age in ages:
    age0=age[0]
    age1=min(age[1],101)
    w=1
    w0=w1=0
    for a in range(age0,age1):
      w0+=w
      w1+=w*dr[a]
      w*=1-dr[a]
    adr[desc][age]=w1/w0

print("      Date      D      E     E/D   D/C*28   E/C*28")
for day in range(datetoday(mindate),datetoday(deaths[-1]['date'])-1):
  caseind=day-datetoday(mcases[0]['date'])
  deathind=day-datetoday(deaths[0]['date'])
  totdeathind=day-datetoday(totdeaths[0]['date'])
  # Estimate background number of deaths that would have occurred in the confirmed cases from day-27 to day
  totc=totd=0
  for (xcases,dr) in [(mcases,adr['m']),(fcases,adr['f'])]:
    for back in range(28):
      for age in ages:
        totc+=xcases[caseind-back][age]
        totd+=xcases[caseind-back][age]*dr[age]/365.2425
  td=totdeaths[totdeathind]['deaths']
  print(daytodate(day),"%6d  %5.1f   %4.1f%%   %5.3f%%   %5.3f%%"%(td,totd,totd/td*100,td/totc*28*100,totd/totc*28*100))
print("      Date      D      E     E/D   D/C*28   E/C*28")
print()
print("D = Number of Covid deaths on this day")
print("E = Expected number of deaths by natural casues on this day from the confirmed cases within the last 28 days")
print("C = Confirmed cases within last 28 days")
print("D/C*28 = CFR")
print("E/C*28 = Natural-cause-IFR")
