from stuff import *
import numpy as np
from math import exp,log,sqrt
from scipy.optimize import minimize
from scipy.special import gammaln,digamma

mindate='2021-05-15'

ladcsv=loadcsv('Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv')
sanger=loadcsv('lineages_by_ltla_and_week.tsv',sep='\t')

lad2reg=dict(zip(ladcsv['LAD19CD'],ladcsv['RGN19NM']))
regions=sorted(list(set(ladcsv['RGN19NM'])))
nreg=len(regions)
maxdate=max(sanger['WeekEndDate'])
ndates=(datetoday(maxdate)-datetoday(mindate))//7+1

# Alpha, Delta maps from region to date to count
A={r:np.zeros(ndates,dtype=int) for r in regions}
D={r:np.zeros(ndates,dtype=int) for r in regions}
for (date,lad,lin,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
  if date>=mindate:
    reg=lad2reg[lad]
    d=ndates-1+(datetoday(date)-datetoday(maxdate))//7
    if lin=='B.1.1.7': A[reg][d]+=n
    if lin=='B.1.617.2': D[reg][d]+=n

def mumax(r,a,d,e,p):
  return gammaln(a+r)+gammaln(d+e*r)-gammaln(r)-gammaln(e*r)+r*(1+e)*log(1-p)

def mumaxd(r,a,d,e,p):
  return digamma(a+r)+e*digamma(d+e*r)-digamma(r)-e*digamma(e*r)+(1+e)*log(1-p)

# q = p/(1-p) = mu/r
def NLL(xx,reg):
  t0,lam,q=xx
  p=q/(1+q)
  LL=0
  for week in range(ndates):
    e=exp(lam*(week*7-t0))
    a,d=A[reg][week],D[reg][week]
    assert a>0 or d>0
    # Need to set r to maximise this:
    # Gamma(a+r)Gamma(d+e*r)/(Gamma(r)*Gamma(e*r))*(1-p)**(r*(1+e))
    r0=-1/((1+e)*log(1-p))
    r1=-1/((1+e)*log(1-p))*(a+d)
    assert mumaxd(r0,a,d,e,p)>=0 and mumaxd(r1,a,d,e,p)<=0
    while r1-r0>1e-12*(r0+r1):
      r=(r0+r1)/2
      if mumaxd(r,a,d,e,p)>=0: r0=r
      else: r1=r
    r=(r0+r1)/2
    LL+=gammaln(a+r)   - gammaln(a+1) - gammaln(r)     + r*log(1-p)+a*log(p)
    LL+=gammaln(d+e*r) - gammaln(d+1) - gammaln(e*r) + e*r*log(1-p)+d*log(p)
  return -LL

def Fisher(xx,reg,eps=1e-4):
  t0,lam,q=xx
  fi=(NLL([t0,lam-eps,q],reg)-2*NLL([t0,lam,q],reg)+NLL([t0,lam+eps,q],reg))/eps**2
  zconf=1.96
  return zconf/sqrt(fi)

for reg in regions:
  res=minimize(NLL,[0,0.1,1],args=(reg,),bounds=[(-50,50), (-0.2,0.2), (1e-6,10)], method="SLSQP")
  if not res.success: raise RuntimeError(res.message)
  t0,lam,q=res.x
  dlam=Fisher(res.x,reg)
  print("%-30s   %6.2f    %5.3f (%5.3f - %5.3f)     %8.3f"%(reg,t0,lam,lam-dlam,lam+dlam,-res.fun),q)
