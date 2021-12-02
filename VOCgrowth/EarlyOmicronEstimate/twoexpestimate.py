from stuff import *
from math import exp,log
from scipy.optimize import minimize
import numpy as np

cases=loadcsv('SAcasecounts.csv')
prov='Total'

N=len(cases['Date'])
day0=datetoday(cases['Date'][0])
day1=datetoday('2021-11-01')
monday=datetoday('2021-11-01')
np.set_printoptions(precision=4,linewidth=1000)

#     Delta            Omicron
# log(exp(x0+x1*(day-day1))+exp(x2+x3*(day-day1))) + x_{4+day%7}
# x1<0, x3>0
# day1 = 2021-11-01 (approximate minimum of cases)
# x0 ~= 5
# x2 ~= 2
# x10  = 1 (gauge fix)

def expand(xx):
  l=[]
  for d in range(N):
    n_delta=exp(xx[0]+xx[1]*(d+day0-day1))
    n_omicron=exp(xx[2]+xx[3]*(d+day0-day1))
    logn=log(n_delta+n_omicron)+xx[4+(day0+d-monday)%7]
    l.append((n_delta,n_omicron,xx[4+(day0+d-monday)%7],logn))
  return l

def NLL(xx,targ):
  l=expand(xx)
  s=0
  for d in range(N):
    s+=(l[d][3]-targ[d])**2/2
  return s

xx=[5, -0.026, 2, 0.22, 0,0,0,0,0,0,0]
bounds=[(xx[0]-5,xx[0]+5), (-0.5,0.1), (xx[2]-5,xx[2]+5), (0.05, 0.4), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (0,0)]
targ=[log(x) for x in cases[prov]]
res=minimize(NLL,xx,args=(targ,),bounds=bounds,method="SLSQP",options={'maxiter':10000})#, 'eps':1e-4, 'ftol':1e-12})

#if not res.success: raise RuntimeError(res.message)
print(res.message)

l=expand(res.x)
for d in range(N):
  print(daytodate(day0+d),"%8.1f  %8.1f  %6.3f  %6.3f   %6.3f"%(l[d][0],l[d][1],l[d][3],targ[d],l[d][3]-targ[d]))
