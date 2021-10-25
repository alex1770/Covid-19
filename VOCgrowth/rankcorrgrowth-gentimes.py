# Trying to find (rho,g) = (growth rate difference, gen time ratio) of two variants using non-parametric methods (Kendall tau - could improve by weighting).
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
  S0=[];S1=[]
  for (dt,v0,v1) in zip(DT,V0,V1):
    sc=scases[dt]/(v0+v1)
    S0.append(sc*v0)
    S1.append(sc*v1)
    print(daytodate(dt),"%8d %8d %8.1f %8.1f"%(v0,v1,sc*v0,sc*v1),file=fp)

# Do simple regression to get good initial values
A=np.array(V0)
D=np.array(V1)
W=1/(1/A+1/D)
day0=int(round(sum(DT)/len(DT)))
X=np.array(DT)-day0
Y=np.log((D+1e-20)/(A+1e-20))
m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
r=np.array([sum(W*Y),sum(W*X*Y)])
c=np.linalg.solve(m,r)
growth=c[1]
day0=day0-c[0]/c[1]

# Not clear whether it's better to scale up to overall samples or not
(N0,N1)=(S0,S1)
#(N0,N1)=(V0,V1)

# Trying to find maximum (over rho) rank correlation of log(n1(t))-rho*log(n0(t)), n0(t), n1(t) = numbers of each variant
# The best rho should be the ratio of gen times, T0/T1
# Rank method using maximum rank correlation doesn't work so well if cases of new variant are just growing, because then it's happy to use rho=0
if 0:
  for i in range(101):
    rho=i/100*2
    l=[log(n1)-rho*log(n0) for (n0,n1) in zip(N0,N1)]
    res=kendalltau(l,range(ndays))
    print("%8.5f  %10.7f"%(rho,res.correlation))

# Trying to find minimum (over rho, g) rank correlation of log(n1(t))-rho*log(n0(t))-g*t, n0(t), n1(t) = numbers of each variant
for g in np.arange(max(growth-.01,0),growth+.01,0.0005):
  for rho in np.arange(0.,2,0.02):
    l=[log(n1)-rho*log(n0)-g*(day-day0) for (day,n0,n1) in zip(DT,N0,N1)]
    res=kendalltau(l,range(ndays))
    if res.pvalue>.025:
      print("%8.5f  %8.5f   %10.7f  %10.7f"%(g,rho,res.correlation,res.pvalue))

