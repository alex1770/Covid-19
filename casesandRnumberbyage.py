import time,calendar,os,json,sys,datetime
from requests import get
from math import sqrt,log,exp
from stuff import *

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_').replace(' to ','_')# cater for 65_to_69, 65_69, "65 to 69" formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

def get_data(req,metric):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  response = get(url+req+'&metric='+metric, timeout=10)
  if not response.ok: raise RuntimeError('Request failed: '+response.text)
  date=time.strftime('%Y-%m-%d',time.strptime(response.headers['Last-Modified'],'%a, %d %b %Y %H:%M:%S %Z'))# Not currently used
  data=response.json()['body']
  
  # Convert from list form to dictionary keyed by age
  day=datetoday(data[0]['date'])
  n=1
  while n<len(data) and datetoday(data[n]['date'])==day-n: n+=1# Find maximal contiguous date range
  data1=[]
  for i in range(n-1,-1,-1):
    d=data[i]
    e={}
    for y in d[metric]:
      if 'cases' in y: val=y['cases']
      else: val=y['deaths']
      try:
        a=parseage(y['age'])
        e[a]=e.get(a,0)+val
      except ValueError:
        pass
    data1.append(e)
  return data1,datetoday(data[n-1]['date'])

req='areaType=nation&areaName=England','newCasesBySpecimenDateAgeDemographics'# England
#req='areaType=utla&areaCode=E08000001','newCasesBySpecimenDateAgeDemographics'# Bolton
#req='areaType=utla&areaCode=E06000008','newCasesBySpecimenDateAgeDemographics'# Blackburn with Darwen
newcases,minday=get_data(*req)
#discard=2# Last entries are unreliable.
#newcases=newcases[:-discard]
n=len(newcases)

def removesupersets(a0):
  a1=[]
  for x in a0:
    ok=1
    for y in a0:
      if y!=x and x[0]<=y[0] and x[1]>=y[1]: ok=0;break
    if ok: a1.append(x)
  a1.sort()
  return a1

ages0=sorted(list(newcases[0].keys()))
baseages=removesupersets(ages0)
ages=[]
for j in range(10):
  ages.append((10*j,10*(j+1)))
for i in range(n):
  for j in range(9):
    newcases[i][(10*j,10*j+10)]=newcases[i][(10*j,10*j+5)]+newcases[i][(10*j+5,10*j+10)]
    newcases[i][(90,100)]=newcases[i][(90,150)]
agestr="       "
for a in ages:
  agestr+="  %4d"%a[0]
agestr+="  %4d"%ages[-1][1]

mgt=5

# One of these two ought to be a multiple of 7
ave=4
gap=7

print(agestr)
for i in range(gap+ave-1,n):
  print(daytodate(minday+i),end='  ')
  for a in ages:
    A=sum(newcases[i-gap-j][a] for j in range(ave))
    B=sum(newcases[i-j][a] for j in range(ave))
    if A>=10: print(" %5.2f"%((B/A)**(mgt/gap)),end='')
    else: print(" -----",end='')
  print()
print(agestr)
print()

print(agestr)
ave=7
Tmax=-1
for i in range(ave-1,n):
  print(daytodate(minday+i),end=' ')
  T=0
  for a in ages:
    B=sum(newcases[i-j][a] for j in range(ave))
    T+=B
    print(" %5.0f"%(B/ave),end='')
  print("   : %6.0f"%(T/ave))
  if daytodate(minday+i)>="2021-05-01" and T>Tmax: Tmax=T;imax=i
print(agestr)
print()

print("Cumulative cases since peak Delta day")
print(agestr)
tot={a:0 for a in ages}
for i in range(imax,n):
  print(daytodate(minday+i),end=' ')
  for a in ages:
    tot[a]+=newcases[i][a]
    print(" %5d"%tot[a],end='')
  print("   : %6d"%(sum(tot.values())))
print(agestr)
tot=0
print(daytodate(minday+n-1),end=' ')
for a in ages:
  t=sum(newcases[i][a] for i in range(n))
  print(" %5d"%t,end='')
  tot+=t
print("   : %6d (cumulative cases over whole pandemic)"%tot)
print(agestr)
print()
