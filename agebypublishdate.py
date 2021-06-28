import time,calendar,os,json,sys,datetime
from stuff import *

# Go back this number of days
nprev=20

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

now=max(os.listdir('apidata'))
casedata=[]
for day in range(datetoday(now)-nprev,datetoday(now)+1):
  dt=daytodate(day)
  with open('apidata/'+dt,'r') as fp: td=json.load(fp)
  d={}
  for x in td[-1]:
    if x=='date': assert td[-1][x]==dt
    else: d[parseage(x)]=td[-1][x]
  casedata.append(d)

ages=sorted(list(casedata[-1]))

for ph in [0,1]:
  print("     Age:",end='')
  for a in ages: print("%3d  "%a[0],end=' ')
  print()
  n=len(casedata)
  for i in range(1,n):
    print(daytodate(datetoday(now)+i-n+1),end=' ')
    t=sum(casedata[i].values())-sum(casedata[i-1].values())
    if ph==0: q=1
    else: q=t/1000
    for a in ages:
      print("%5d"%((casedata[i][a]-casedata[i-1][a])/q+.5),end=' ')
    print("= %d"%(t/q+.5))
  print()
