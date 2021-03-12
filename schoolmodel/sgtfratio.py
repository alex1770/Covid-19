import numpy as np
import csv
from math import log

fn="S-cases-age-stp-week-nhser-labeled-1446-PHE.csv"

np.set_printoptions(precision=6,linewidth=120)

d={}
e={}
weeks=set()
nhs=set()
locs=set()
with open(fn,"r") as fp:
  r=csv.reader(fp)
  headings=next(r)
  for row in r:
    week=int(row[0]);weeks.add(week)
    loc=row[1];locs.add(loc)
    nhser=row[2];nhs.add(nhser)
    ageband=row[3]
    agenum=int(row[3].split('-')[0])//10
    neg=float(row[4])
    pos=float(row[5])
    tot=int(row[9])
    d[week,loc,agenum]=(neg,pos)
    if (week,agenum) not in e: e[week,agenum]=[0,0]
    e[week,agenum][0]+=neg
    e[week,agenum][1]+=pos

weeks=sorted(list(weeks))
locs=sorted(list(locs))
nhs=sorted(list(nhs))

for t in range(1,8):
  w0,w1=weeks[0],weeks[-1]
  group0=range(0,t)
  group1=range(t,8)
  acc=k=0
  for loc in locs:
    if 0:
      ok=[(w,loc,a) in d for w in [w0,w1] for a in range(8)]
      if not all(ok): continue
      r0=sum(d[w0,loc,a][0] for a in group0)/sum(d[w0,loc,a][1] for a in group0)
      r1=sum(d[w0,loc,a][0] for a in group1)/sum(d[w0,loc,a][1] for a in group1)
      r2=sum(d[w1,loc,a][0] for a in group0)/sum(d[w1,loc,a][1] for a in group0)
      r3=sum(d[w1,loc,a][0] for a in group1)/sum(d[w1,loc,a][1] for a in group1)
      if r0*r3>0 or r1*r2>0:
        acc+=1
        # r=r0/r1/(r2/r3)
        # print(log(r))
        if r0*r3<=r1*r2: k+=1
    if 1:
      ok=[(w,loc,a) in d for w in [w1] for a in range(8)]
      if not all(ok): continue
      r2=sum(d[w1,loc,a][0] for a in group0)/sum(d[w1,loc,a][1] for a in group0)
      r3=sum(d[w1,loc,a][0] for a in group1)/sum(d[w1,loc,a][1] for a in group1)
      acc+=1
      if r2<=r3: k+=1
    if 0:
      ok=[(w,loc,a) in d for w in weeks for a in range(8)]
      if not all(ok): continue
      r2=sum(d[w,loc,a][0] for a in group0 for w in weeks)/sum(d[w,loc,a][1] for a in group0 for w in weeks)
      r3=sum(d[w,loc,a][0] for a in group1 for w in weeks)/sum(d[w,loc,a][1] for a in group1 for w in weeks)
      acc+=1
      if r2<=r3: k+=1

  rej=len(locs)-acc
  print()
  print("Separating into age groups %d-%d and %d-%d"%(0,t*10,t*10,80))
  print("Locations accepted:",acc)
  print("Locations rejected:",rej)
  print("Number of locations where lower age group has more SGTF:",acc-k,"out of",acc)
  totp=0;p=1;n=acc
  for r in range(k+1):
    totp+=p
    p*=(n-r)/(r+1)
  print("p = %.3g"%(totp/2**n))
  
