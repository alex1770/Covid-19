# Starting from ONS dropout data, try to make model predicting proportions of dropouts of OR, N, S.
# Following on from Theo Sanderson's analysis at https://theo.io/post/2021-01-22-ons-data/.
# Data from tabs 6a, 6b of spreadsheet here https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveypilot/22january2021/relateddata
# as transcribed at https://github.com/theosanderson/theo.io/blob/master/content/post/2021-01-22-ons-data/ons_ct.csv
# and updated with data from subsequent ONS infection surveys.
# This differs from dropoutmodel.py in that it doesn't assume constant logistic growth.

# Model:
# Let logistic(x)=exp(x)/(1+exp(x))
# Relative prevalence of B.1.1.7 = logistic(logodds(r,t)), where r=region and logodds(r,t) has approx constant growth in t (see below for how logodds(r,t) is defined in terms of parameters)
# Choose a uniform random number Z in [-5,5], which is fixed for the three genes, representing viral load,
# then the probability of dropout for gene X (N, OR or S) = logistic(b*(Ct-Z-a_X)).
# The probability of SGTF (mandatory S dropout, assumed to be due to B.1.1.7) is exp(logodds[r,t]) with logodds[r,t] given as below:
# Parameters (4+2*nregions+ndates-2):
#   a_X        : 3 parameters, one for each gene N, OR and S, encoding their "robustness" (lower = more fragile)
#   b          : 1 parameter ("ctmult") encoding dependence of dropout probability on Ct
#   logodds0   : nregions parameters encoding logodds in region r at time index 0
#   dlogodds0  : nregions parameters encoding dlogodds in region r at time index 0
#   ddlogodds  : ndates-2 parameters encoding difference in dlogodds at time index t (independent of region)
# which expand to dlogodds[r,t] and logodds[r,t] (t=0,...,ndates-1; r=0,...,nregions-1) satisfying:
#   dlogodds[r,0] = dlogodds0[r]
#   dlogodds[r,t+1] - dlogodds[r,t] = ddlogodds[t] (t=0,...,ndates-2)
#   logodds[r,0] = logodds0[r]
#   logodds[r,t+1]-logodds[r,t] = dlogodds[r,t] (t=0,...,ndates-2)
#
# Idea is that the 2nd difference (ddlogodds) would be zero if simple constant logistic growth applied, so treat 2nd diff as a perturbation.
# Expect logistic growth to slow a little under lockdown conditions. These are (currently) England-wide so 2nd diffs are region-independent,
# which also keeps the free parameter count down.
# (Allowing ddlogodds to depend on r was tried and doesn't improve the fit a huge amount.)

import sys,time,calendar,csv
from math import log,exp,sqrt
import numpy as np
from scipy.optimize import minimize

mindate="2020-10-01"

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

class testdata:
  def __init__(self,date):
    self.date=date
    self.t=datetoday(date)-datetoday(date0)
    self.p=np.zeros([2,2,2])# p[whether N][whether OR][whether S] = proportion (adding to 1)
    self.Ct=0
  def __repr__(self):
    s=self.date+" :"
    for t in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(1,1,1)]: s+=" %5.1f"%(self.p[t]*100)
    s+=" : %.1f"%(self.Ct)
    return s

# Convert flat array into structured parameters (lvalues)
def paramnames(xx):
  robustness=xx[:3]
  ctmult=xx[3:4]
  logodds0=xx[4:4+nregions]
  dlogodds0=xx[4+nregions:4+2*nregions]
  ddlogodds=xx[4+2*nregions:]# length ndates-2
  return (robustness,ctmult,logodds0,dlogodds0,ddlogodds)

def getlogodds(xx):
  (robustness,ctmult,logodds0,dlogodds0,ddlogodds)=paramnames(xx)
  logodds=np.zeros((nregions,ndates),dtype=float)
  logodds[:,0]=logodds0
  dlogodds=dlogodds0.copy()
  for t in range(ndates-1):
    logodds[:,t+1]=logodds[:,t]+dlogodds
    if t<ndates-2: dlogodds+=ddlogodds[t]
  return logodds

# csv headings:
# 0             1       2               3       4       5       6       7       8       9       10
# RegionType	Region	Week started	N only	OR only	S only	OR+N	OR+S	N+S	OR+N+S	Mean
#
# 11               12                   13              14              15
# 10th Percentile  25th Percentile	50th Percentile	75th Percentile	90th Percentile

data={}
date0=None
with open("ons_ct.csv","r") as fp:
  reader=csv.reader(fp)
  headings=next(reader)
  for row in reader:
    if row[0]=="EnglandRegion":
      date=time.strftime("%Y-%m-%d",time.strptime(row[2],"%d %B %Y"))
      if date>=mindate:
        if date0==None: date0=date
        d=testdata(date)
        assert d.t>=0 and d.t%7==0
        d.p[1][0][0]=float(row[3])# N
        d.p[0][1][0]=float(row[4])# OR
        d.p[0][0][1]=float(row[5])# S
        d.p[1][1][0]=float(row[6])# OR+N
        d.p[0][1][1]=float(row[7])# OR+S
        d.p[1][0][1]=float(row[8])# N+S
        d.p[1][1][1]=float(row[9])# OR+N+S
        d.p=d.p/d.p.sum()
        d.Ct=float(row[10])# Mean Ct
        data.setdefault(row[1],[]).append(d)

regions=sorted(list(data))
nregions=len(regions)
x=set(len(x) for x in data.values());assert len(x)==1
ndates=x.pop()
nparams=3+1+nregions*2+ndates-2
smoothness=5.0

# Work out expected dropout matrix for region r, dropout matrix, date and Ct value d, model parameters xx[]
def estimatedropoutmatrix(r,d,robustness,ctmult,logodds):
  tc=np.zeros([2,2,2])
  for offset in range(-5,6):# Integrate over viral load
    p=1/(1+exp(-logodds[r][d.t//7]))# Relative prevalence of B.1.1.7
    dp=[1/(1+exp(-ctmult[0]*(d.Ct-a+offset*1.0))) for a in robustness]# Probability of dropout for each gene
    dp[2]=p+(1-p)*dp[2]# Can treat B.1.1.7 as something that increases probability of S gene dropout
    c=np.zeros([2,2,2])
    for i in range(2):
      for j in range(2):
        for k in range(2):
          c[i,j,k]=(dp[0] if i==0 else 1-dp[0])*(dp[1] if j==0 else 1-dp[1])*(dp[2] if k==0 else 1-dp[2])
    tc+=c
  tc[0,0,0]=0;tc=tc/tc.sum()# Remove all-dropout and renormalise
  return tc

# Calculate modelling error (cross entropy) associated with parameters xx[]
# Would really like to multiply cross entropy for each (region,week) by the number of tests done, but that information is not available.
def err(xx):
  (robustness,ctmult,logodds0,dlogodds0,ddlogodds)=paramnames(xx)
  logodds=getlogodds(xx)
  E=0
  for (r,region) in enumerate(regions):
    for d in data[region]:
      c=estimatedropoutmatrix(r,d,robustness,ctmult,logodds)
      # -log(likelihood of actual dropout matrix if true probs are from estimated dropout matrix)
      E-=(d.p*np.log(c+1e-100)).sum()
  return E+smoothness*(ddlogodds*ddlogodds).sum()

# Baseline error to convert cross entropy into KL divergence
def err0():
  E=0
  for (r,region) in enumerate(regions):
    for d in data[region]:
      E-=(d.p*np.log(d.p+1e-100)).sum()
  return E

# Initial parameter values and bounds
xx=np.zeros(nparams);(robustness,ctmult,logodds0,dlogodds0,ddlogodds)=paramnames(xx)
lbound=np.zeros(nparams);(lbound_r,lbound_c,lbound_l,lbound_d0,lbound_dd)=paramnames(lbound)
ubound=np.zeros(nparams);(ubound_r,ubound_c,ubound_l,ubound_d0,ubound_dd)=paramnames(ubound)
robustness.fill(28);lbound_r.fill(10);ubound_r.fill(50)
ctmult[0]=1.5;lbound_c[0]=-1;ubound_c[0]=3
logodds0.fill(-10);lbound_l.fill(-30);ubound_l.fill(30)
dlogodds0.fill(1.0);lbound_d0.fill(-1);ubound_d0.fill(1)
ddlogodds.fill(0.0);lbound_dd.fill(-1);ubound_dd.fill(1)
bounds=list(zip(lbound,ubound))

print("Initial total KL divergence + prior on 2nd diffs = %.1f bits"%((err(xx)-err0())/log(2)))
print("Using smoothness coefficient %.3f"%smoothness)

res=minimize(err,xx,method="SLSQP",bounds=bounds,options={"maxiter":1000})
print(res.message)
if not res.success: sys.exit(1)

xx=res.x
(robustness,ctmult,logodds0,dlogodds0,ddlogodds)=paramnames(xx)
logodds=getlogodds(xx)
KL=res.fun-err0()
n=sum(len(data[region]) for region in regions)
print("Total KL divergence + prior on 2nd diffs = %.1f bits (this is missing factors of number of tests done because that information isn't published, so understates the information deficit)"%(KL/log(2)))
print("Average KL divergence + prior on 2nd diffs = %.3f bits (ditto)"%(KL/n/log(2)))
print()

now=time.strftime('%Y-%m-%d',time.localtime())
day0=datetoday(date0)
tnow=datetoday(now)-day0

# Determine max growth
gmax=np.max(logodds[:,1:]-logodds[:,:-1],axis=1)/7

# Interpolate/extrapolate logodds
logodds_i=np.zeros((nregions,tnow+1),dtype=float)
for r in range(nregions):
  g=(logodds[r][ndates-1]-logodds[r][ndates-2])/(1*7)
  for t in range(tnow+1):
    if t<7*(ndates-1):
      l=((7-t%7)*logodds[r][t//7]+(t%7)*logodds[r][t//7+1])/7
    else:
      l=logodds[r][ndates-1]+(t-7*(ndates-1))*g
    logodds_i[r][t]=l

print("Robustness of N  = %.1f"%robustness[0])
print("Robustness of OR = %.1f"%robustness[1])
print("Robustness of S  = %.1f"%robustness[2])
print("Dependence of dropout on Ct = %.3f"%ctmult[0])
print()
print("Region                    Est'd crossover     Extrapolated %relative")
print("------                    date of B.1.1.7     prevalence on",now,"   Maximum growth rate")
print("                          ---------------     ------------------------    -------------------")
for (r,region) in enumerate(regions):
  for t in range(tnow):
    if logodds_i[r][t]<0 and logodds_i[r][t+1]>=0: tc=t+logodds_i[r][t+1]/(logodds_i[r][t+1]-logodds_i[r][t]);break
  else: tc=None
  p=1/(1+exp(-logodds_i[r][tnow]))
  print("%-24s  %s        %6.1f                       %6.3f"%(region,daytodate(day0+tc) if tc!=None else "????-??-??",p*100,gmax[r]))

print()

with open('dropoutmodel.csv','w') as fp:
  print(", ".join(["date"]+regions),file=fp)
  for t in range(tnow+1):
    print(daytodate(day0+t),end="",file=fp)
    for r in range(len(regions)):
      p=1/(1+exp(-logodds_i[r][t]))
      print(", %10.4g"%(p),end="",file=fp)
    print(file=fp)

def printdropouts(m):
  for t in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(1,1,1)]:
    print(" %5.1f"%(m[t]*100),end="")

print("                                                           Actual                                        Estimate                 ")
print("                                         ------------------------------------------     ------------------------------------------")
print("                  Region       Date         N    OR     S  OR+N  OR+S   N+S  OR+N+S        N    OR     S  OR+N  OR+S   N+S  OR+N+S")
(robustness,ctmult,logodds0,dlogodds0,ddlogodds)=paramnames(xx)
logodds=getlogodds(xx)
for (r,region) in enumerate(regions):
  for d in data[region]:
    print("%24s"%region,d.date,end="    ")
    printdropouts(d.p)
    print("     ",end="")
    c=estimatedropoutmatrix(r,d,robustness,ctmult,logodds)
    printdropouts(c)
    print()
