import json
import numpy as np
from stuff import *

with open('SouthAfricaHospData.json') as fp:
  hd=json.load(fp)

provinces=['Eastern Cape', 'Free State', 'Gauteng', 'KwaZulu-Natal', 'Limpopo', 'Mpumalanga', 'North West', 'Northern Cape', 'Western Cape']

mindate='2021-02-03'
minday=datetoday(mindate)
maxday=datetoday(max(hd))
n=maxday-minday+1

for day in range(minday-1,maxday):
  date0=daytodate(day)
  date1=daytodate(day+1)
  print(date0,end='')
  for key in 'Admissions to Date', 'Died to Date':
    print('  ',end='')
    for prov in provinces:
      d=hd[date1][prov][key]-hd[date0][prov][key]
      if d<0: print('X',end='')
      else: print('.',end='')
    prov='South Africa'
    d=hd[date1][prov][key]-hd[date0][prov][key]
    if d<0: print('*',end='')
    else: print(' ',end='')
  print()
  
aa=[];bb=[]
loc='South Africa'
for day in range(minday-1,maxday+1):
  date=daytodate(day)
  if date not in hd: print(date,"missing");continue
  aa.append(hd[date][loc]['Admissions to Date'])
  bb.append(hd[date][loc]['Died to Date'])
  #print(date,aa[-1],bb[-1])

for i in range(n):
  if aa[i+1]<aa[i]: print(daytodate(minday-1+i),"aa anomaly",aa[i+1]-aa[i],aa[i]-aa[i-1])
  if bb[i+1]<bb[i]: print(daytodate(minday-1+i),"bb anomaly",bb[i+1]-bb[i],aa[i]-aa[i-1])
  
for d in range(datetoday('2021-11-15'),datetoday('2021-11-26')):
  date=daytodate(d)
  print(date,end='')
  for key in 'Admissions to Date', 'Died to Date':
    print('       ',end='')
    for prov in provinces:
      print(" %6d"%hd[date][prov][key],end='')
  print()
print()

for d in range(datetoday('2021-11-15'),datetoday('2021-11-26')):
  date=daytodate(d)
  date0=daytodate(d-1)
  date1=daytodate(d+1)
  print(date,end='')
  for key in 'Admissions to Date', 'Died to Date':
    print('       ',end='')
    for prov in provinces:
      x=hd[date][prov][key]-(hd[date0][prov][key]+hd[date1][prov][key])/2
      y=(hd[date0][prov][key]-hd[daytodate(d-21)][prov][key])/20
      if y!=0: print(" %6.1f"%(x/y),end='')
      else: print("      -",end='')
  print()
print()
