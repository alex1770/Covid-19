# Trying to find (rho,g) = (growth rate difference, gen time ratio) of two variants using simple weighted regression
# It relies on non-steady overall growth in both variants over the time window to get a localised result in (rho,g)-space.
# Steady overall growth in one variant means you can trade off increased coefficient of that variant with increased g.

from stuff import *
import numpy as np
from math import exp,log,sqrt
from scipy.stats import kendalltau

mindate='0000-00-00'
maxdate='9999-99-99'
mincount=10
prlevel=1
if len(sys.argv)>1: mindate=sys.argv[1]
if len(sys.argv)>2: maxdate=sys.argv[2]
if len(sys.argv)>3: mincount=int(sys.argv[3])
print("Using date range",mindate,"-",maxdate,file=sys.stderr)
print("Min count:",mincount,file=sys.stderr)

# Variant0, Variant1 counts by day
V0=[];V1=[];DT=[]
fp=sys.stdin
#fp=open('cog.y145h','r')
#fp=open('alphadelta','r')
if 1:
  for x in fp:
    y=x.strip().split()
    if y[0]>=mindate and y[0]<=maxdate:
      d=datetoday(y[0])
      v0=int(y[1])
      v1=int(y[2])
      if v0>=mincount and v1>=mincount: V0.append(v0);V1.append(v1);DT.append(d)
ndays=len(V0)

# Scale sequenced totals up to smoothed version of overall samples (attempting to correct for limited sequencing capacity and testing fluctuations)
casescsv=loadcsv('UKcasesbysampledate.csv')
cases=dict(zip( map(datetoday,casescsv['date']) , casescsv['newCasesBySpecimenDate'] ))
scases={}
for day in range(min(cases),max(cases)+1):
  s0=s1=0
  for day1 in range(day-3,day+4):
    if day1 in cases: s0+=1;s1+=cases[day1]
  if s0>0: scases[day]=s1/s0

with open('temp1','w') as fp:
  S0=[];S1=[];SC=[]
  for (dt,v0,v1) in zip(DT,V0,V1):
    sc=scases[dt]/(v0+v1)
    S0.append(sc*v0)
    S1.append(sc*v1)
    SC.append(sc)
    print(daytodate(dt),"%8d %8d %8d %8.1f %8.1f"%(dt,v0,v1,sc*v0,sc*v1),file=fp)

# Not clear whether it's better to scale up to overall samples or not
(N0,N1,SF)=(S0,S1,np.array(SC))
#(N0,N1,SF)=(V0,V1,np.zeros(ndays)+1)

N0=np.array(N0)
N1=np.array(N1)

zconf=1.96
for rho in np.arange(0,2,0.1):#[0,1,10,1e4]:#np.arange(0,10,1):
  W=1/((rho**2/N0+1/N1)*SF*SF)# inverse variance
  day0=sum(DT)/len(DT)# Initial guess to keep X nicely conditioned
  X=np.array(DT)-day0
  Y=np.log((N1+1e-20))-rho*np.log((N0+1e-20))
  m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
  r=np.array([sum(W*Y),sum(W*X*Y)])
  c=np.linalg.solve(m,r)
  growth=c[1]
  day0=day0-c[0]/c[1]
  mi=np.linalg.pinv(m)
  R=c[0]+c[1]*X-Y
  resid=(R*R*W).sum()/len(R)
  err=zconf*sqrt(mi[1,1])
  print("%8.5f  %s %8.2f   %5.2f%% (%5.2f%% - %5.2f%%)   %8.2f"%(rho,daytodate(day0+.5),day0,c[1]*100,(c[1]-err)*100,(c[1]+err)*100,resid))

