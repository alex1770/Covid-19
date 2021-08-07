from stuff import *
import numpy as np
from math import exp,log,sqrt
from scipy.optimize import minimize

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

def NLL(xx,reg):
  t0,lam=xx
  LL=0
  for d in range(ndates):
    e=exp(lam*(d*7-t0))
    LL+=A[reg][d]*log(1/(1+e))+D[reg][d]*log(e/(1+e))
  return -LL

def Fisher(xx,reg,eps=1e-4):
  t0,lam=xx
  fi=(NLL([t0,lam-eps],reg)-2*NLL([t0,lam],reg)+NLL([t0,lam+eps],reg))/eps**2
  zconf=1.96
  return zconf/sqrt(fi)

for reg in regions:
  res=minimize(NLL,[0,0.1],args=(reg,),bounds=[(-50,50), (-0.2,0.2)], method="SLSQP")
  if not res.success: raise RuntimeError(res.message)
  t0,lam=res.x
  dlam=Fisher(res.x,reg)
  print("%-30s   %6.2f    %5.3f (%5.3f - %5.3f)"%(reg,t0,lam,lam-dlam,lam+dlam))
