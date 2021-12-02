from stuff import *
from math import exp,log
from scipy.optimize import minimize
import numpy as np

cases=loadcsv('SAcasecounts.csv')
prov='Total'

N=len(cases['Date'])
day0=datetoday(cases['Date'][0])
monday=datetoday('2021-11-01')
np.set_printoptions(precision=4,linewidth=1000)

#     Delta            Omicron
# log(exp(x1*(day-day1)+x0)+exp(x3*(day-x2))) + x_{4+day%7}
# x1<0, x3>0
# x0 ~~= 2022-03-08
# x2  ~= 2021-10-25
# x10  = 1 (gauge fix)

def NLL(xx,targ,deb=False):
  s=0
  for d in range(N):
    n_delta=exp(xx[1]*(d-xx[0]))
    n_omicron=exp(xx[3]*(d-xx[2]))
    l=log(n_delta+n_omicron)+xx[4+(day0+d-monday)%7]
    if deb:
      print(daytodate(day0+d),n_delta,n_omicron,l,targ[d],l-targ[d])
    s+=(l-targ[d])**2/2
  return s

xx=[datetoday('2022-03-01')-day0, -0.026, datetoday('2021-10-20')-day0, 0.22, 0,0,0,0,0,0,0]
bounds=[(xx[0]-30,xx[0]+30), (-0.5,0.1), (xx[2]-30,xx[2]+30), (0.05, 0.4), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (0,0)]
targ=[log(x) for x in cases[prov]]
res=minimize(NLL,xx,args=(targ,),bounds=bounds,method="SLSQP",options={'maxiter':10000})#, 'eps':1e-4, 'ftol':1e-12})

#if not res.success: raise RuntimeError(res.message)
print(res.message)

NLL(res.x,targ,True)
print()

