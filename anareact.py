# Analysing growth trend from REACT-1 survey, published on 21 Jan 2021, mostly covering period 5-15 Jan 2021
# Data from https://github.com/mrc-ide/reactidd

import os,time,calendar,csv
from math import log,exp,sqrt
from scipy.optimize import minimize

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def load(fn,regions,mindate="0000-00-00"):
  with open(os.path.join(ddir,fn),'r') as fp:
    reader=csv.reader(fp)
    locs=next(reader)[1:]
    l=[]
    for row in reader:
      if row[0]>=mindate:
        l.append((datetoday(row[0]),sum(int(row[i]) for i in regions)))
  return l

ddir="reactidd/inst/extdata"

# -log likelihood, with a growth parameter
# xx=[p/(1-p), growth]
def LL(xx,pos,tot):
  (lam,g)=xx# lam=p/(1-p), the odds
  n=len(pos)
  lam*=exp(-g*(n-1)/2)# Centre so no growth factor in the middle (WLOG but may help optimisation)
  glam=exp(g)
  s=0
  for (np,nt) in zip(pos,tot):
    s+=np*log(lam)-nt*log(1+lam)
    lam*=glam
  return -s

# Work out Fisher information with respect to g
def Fisher(lam,g,pos,tot):
  eps=1e-3
  L0=LL([lam,g-eps],pos,tot)
  L1=LL([lam,g    ],pos,tot)
  L2=LL([lam,g+eps],pos,tot)
  return (L0-2*L1+L2)/eps**2

# Regions from csv file:
#    0         1           2           3                   4                  5               6             7            8         9
# [Date,] South East, North East, North West, Yorkshire and The Humber, East Midlands, West Midlands, East of England, London, South West

for (regions,desc) in [(range(1,10),"England"), ([1,7,8,9],"South England"), ([2,3,4,5,6],"North and Central England")]:

  print("REACT-1 survey:",desc)
  print()
  pos0=load("positive.csv",regions,"2020-12-25")
  tot0=load("total.csv",regions,"2020-12-25")
  
  # Match dates in pos and tot and convert to simple list
  minday=min(x for (x,y) in pos0+tot0)
  maxday=max(x for (x,y) in pos0+tot0)
  n=maxday-minday+1
  pos=[0]*n;tot=[0]*n
  for (x,y) in pos0: pos[x-minday]=y
  for (x,y) in tot0: tot[x-minday]=y
  
  for d in range(n):
    print(daytodate(minday+d),"%6d  %6d"%(pos[d],tot[d]))
  print()
  
  # Try no growth
  p=sum(pos)/sum(tot)
  L=sum(pos)*log(p)+(sum(tot)-sum(pos))*log(1-p)
  print("Best constant prevalence (no growth): p=%.2f%%"%(p*100))
  print()

  # Try a growth parameter
  res=minimize(LL,[0.01,0],args=(pos,tot),method="SLSQP",bounds=[(1e-4,10),(-10,10)])
  if not res.success: raise RuntimeError(res.message)
  print("Extra parameter for growth increases log likelihood by %.3f,"%(-res.fun-L))
  (lam,g)=res.x
  p=lam/(1+lam)
  gerr=1/sqrt(Fisher(lam,g,pos,tot))
  print("using p=%.2f%% on %s and growth of %.1f%%/day"%(p*100,daytodate(minday+(n-1)/2),g*100))
  print("2/sqrt(Fisher information of g) = %.1f%%/day (some sort of error estimate for g)"%(2*gerr*100))
  print("\n-----------------------------------------------\n")
