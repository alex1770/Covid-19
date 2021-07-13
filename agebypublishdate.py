import time,calendar,os,json,sys,datetime
from stuff import *

# Go back this number of days
nprev=50

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

newcases=[]
dates=[]
n=len(casedata)
for i in range(1,n):
  dates.append(daytodate(datetoday(now)+i-n+1))
  newcases.append({})
  for a in ages:
    newcases[-1][a]=casedata[i][a]-casedata[i-1][a]

agestr="     Age:"+'   '.join("%3d"%a[0] for a in ages)

print(agestr)
for i in range(n-1):
  print(dates[i],end=' ')
  t=sum(newcases[i].values())
  for a in ages:
    print("%5d"%(newcases[i][a]+.5),end=' ')
  print("= %d"%(t+.5))
print()

print(agestr)
for i in range(n-1):
  print(dates[i],end=' ')
  t=sum(newcases[i].values())
  for a in ages:
    print("%5d"%(newcases[i][a]/t*1000+.5),end=' ')
  print("= 1000")
print()

mgt=5
print(agestr)
ave=7
for i in range(7+ave-1,n-1):
  print(dates[i],end=' ')
  tA=tB=0
  for a in ages:
    A=sum(newcases[i-j][a] for j in range(ave))
    B=sum(newcases[i-7-j][a] for j in range(ave))
    tA+=A;tB+=B
    if B>0:
      print("%5.2f"%((A/B)**(mgt/7)),end=' ')
    else:
      print(" ????",end=' ')
  print(": %5.2f"%(tA/tB))
print(agestr)
print()

extrapr=1
mgt=5
ages1=[(i*10,(i+1+5*(i==9))*10) for i in range(10)]
for c in newcases:
  for (a,b) in ages1:
    c[(a,b)]=sum(c[(a1,b1)] for (a1,b1) in ages if a1>=a and b1<=b)
agestr="     Age:"+'   '.join("%3d"%a[0] for a in ages1)
print(agestr)
ave=3
for i in range(7+ave-1,n-1):
  if extrapr:
    print(dates[i],end=' ')
    tA=tB=0
    for a in ages1:
      A=sum(newcases[i-j][a] for j in range(ave))
      B=sum(newcases[i-7-j][a] for j in range(ave))
      tA+=A;tB+=B
      print("%5d"%A,end=' ')
    print(": %6d"%tA)
    
    print(dates[i],end=' ')
    tA=tB=0
    for a in ages1:
      A=sum(newcases[i-j][a] for j in range(ave))
      B=sum(newcases[i-7-j][a] for j in range(ave))
      tA+=A;tB+=B
      print("%5d"%B,end=' ')
    print(": %6d"%tB)

  print(dates[i],end=' ')
  tA=tB=0
  s=0
  for a in ages1:
    A=sum(newcases[i-j][a] for j in range(ave))
    B=sum(newcases[i-7-j][a] for j in range(ave))
    tA+=A;tB+=B
    if a[0]<=80: s+=(a[0]-40)*A/B
    if B>0:
      print("%5.2f"%((A/B)**(mgt/7)),end=' ')
    else:
      print(" ????",end=' ')
  print(": %5.2f"%(tA/tB),s)

  if extrapr: print()
  
print(agestr)
print()
