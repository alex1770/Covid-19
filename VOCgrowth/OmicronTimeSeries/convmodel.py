import json
import numpy as np
from stuff import *

with open('SouthAfricaHospData.json') as fp:
  data=json.load(fp)

provinces=['Eastern Cape', 'Free State', 'Gauteng', 'KwaZulu-Natal', 'Limpopo', 'Mpumalanga', 'North West', 'Northern Cape', 'Western Cape']

day0=Date(min(data))
now=Date(max(data))
locations=list(data[str(now)])
headings=list(data[str(now)]['South Africa'])
monotonic={'Admissions to Date', 'Died to Date', 'Discharged to Date'}

if 0:
  loc='Gauteng'
  #key='Admissions to Date'
  key='Died to Date'
  for day in Daterange(day0+1,now+1):
    try:
      print(day,data[str(day)][loc][key]-data[str(day-1)][loc][key],data[str(day)][loc][headings[-1]])
    except:
      pass

if 0:
  loc='Gauteng'
  #head='Admissions in Previous Day'
  head='Died to Date'
  for day in Daterange(now,day0,-1):
    if not (str(day) in data and head in data[str(day)][loc]): break
  l=[data[str(d)][loc][head]-data[str(d-1)][loc][head] for d in Daterange(day+2,now+1)]
  m=weekdayadj([max(x,0) for x in l],10)
  for (x,y,z) in zip(Daterange(day+1,now+1),l,m):print(x,y,z)
  #for d in range(7):
  #  print(sum(l[d::7]),sum(m[d::7]))

n=now-day0+1
for loc in locations:
  for head in headings:
    # A day is classified as 'bad' if it doesn't exist or, for a monotonic quantity, if the value is less than something in the past or greater than something in the future
    bad=[0]*n
    for day in Daterange(day0,now+1):
      if str(day) not in data or head not in data[str(day)][loc]: bad[day-day0]=1
    loc='Gauteng';head='Died to Date'#alter
    if head in monotonic:
      mx=-1e9
      for day in Daterange(day0,now+1):
        if str(day) in data and head in data[str(day)][loc]:
          v=data[str(day)][loc][head]
          if v<mx: bad[day-day0]=1
          else: mx=v
      mi=1e9
      for day in Daterange(now,day0-1,-1):
        if str(day) in data and head in data[str(day)][loc]:
          v=data[str(day)][loc][head]
          if v>mi: bad[day-day0]=1
          else: mi=v
    for day in Daterange(day0,now+1):
      if str(day) in data and head in data[str(day)][loc]:
        v=data[str(day)][loc][head]
      else: v='-'
      print(day,bad[day-day0],"%8s"%str(v))
    poi
    #day=day0
    #while day<=now and str(day) not in data: day+=1
    #while day<=now:
    #  # Look for stretches of 'bad' date transitions (involving missing date, or decrease in monotonic quantity)
    #  day1=day
    #  while day<=now and (str(day+1) not in data or (head in monotonic and data[str(day+1)][loc][head]<data[str(day1)][loc][head])): day+=1
    #  if day1<day: print(loc,'/',head,':',str(day1),"-",str(day))
    #  while day<=now and (str(day+1) in data and (head not in monotonic or data[str(day+1)][loc][head]>=data[str(day)][loc][head])): day+=1

oi

mindate='2021-02-03'
minday=Date(mindate)
maxday=Date(max(data))
n=maxday-minday+1

for day in range(minday-1,maxday):
  date0=str(day)
  date1=str(day+1)
  print(date0,end='')
  for key in 'Admissions to Date', 'Died to Date':
    print('  ',end='')
    for prov in provinces:
      d=data[date1][prov][key]-data[date0][prov][key]
      if d<0: print('X',end='')
      else: print('.',end='')
    prov='South Africa'
    d=data[date1][prov][key]-data[date0][prov][key]
    if d<0: print('*',end='')
    else: print(' ',end='')
  print()
  
aa=[];bb=[]
loc='South Africa'
for day in Daterange(minday-1,maxday+1):
  date=str(day)
  if date not in data: print(date,"missing");continue
  aa.append(data[date][loc]['Admissions to Date'])
  bb.append(data[date][loc]['Died to Date'])
  #print(date,aa[-1],bb[-1])

for i in range(n):
  if aa[i+1]<aa[i]: print(str(minday-1+i),"aa anomaly",aa[i+1]-aa[i],aa[i]-aa[i-1])
  if bb[i+1]<bb[i]: print(str(minday-1+i),"bb anomaly",bb[i+1]-bb[i],aa[i]-aa[i-1])

for d in Daterange(Date('2021-11-15'),Date('2021-11-26')):
  date=str(d)
  print(date,end='')
  for key in 'Admissions to Date', 'Died to Date':
    print('       ',end='')
    for prov in provinces:
      print(" %6d"%data[date][prov][key],end='')
  print()
print()

for d in Daterange(Date('2021-11-15'),Date('2021-11-26')):
  date=str(d)
  date0=str(d-1)
  date1=str(d+1)
  print(date,end='')
  for key in 'Admissions to Date', 'Died to Date':
    print('       ',end='')
    for prov in provinces:
      x=data[date][prov][key]-(data[date0][prov][key]+data[date1][prov][key])/2
      y=(data[date0][prov][key]-data[daytodate(d-21)][prov][key])/20
      if y!=0: print(" %6.1f"%(x/y),end='')
      else: print("      -",end='')
  print()
print()
