import requests,sys
import numpy as np
from math import log,exp,sqrt
from stuff import *
from scipy.optimize import minimize
from random import normalvariate as normal,random,seed,randrange
from parseonsprevinc import getonsprevinc

startdate=Date("2021-06-01")
if len(sys.argv)>1: startdate=Date(sys.argv[1])
print("Start date =",startdate)

ignore={"2020-12-25","2020-12-26","2020-12-27","2020-12-28",
        "2021-12-25","2021-12-26","2021-12-27","2021-12-28",
        "2022-12-25","2022-12-26","2022-12-27"}

#ignore={"2020-12-25","2020-12-28","2021-12-25","2021-12-26","2021-12-27","2021-12-28","2022-12-25","2022-12-26","2022-12-27"}
#ignore={}
  
# population=55.98e6
monday=Date("2022-01-03")
scale=1e5# Units of this number of people

np.set_printoptions(precision=3,suppress=True,linewidth=200)

minback=1
maxback=10
odp=4
# Order=1 if you think the prior for xx[] (CAR) is Brownian motion-like (in particular, Markov)
# Order=2 if you think the prior for xx[] is more like the integral of Brownian motion.
# I think xx[] does actually have "momentum" sometimes - there is a downwards (or upwards) trend that is predictive, so perhaps order=2 is better?
# Though it can dart around sometimes for random reasons, like Christmas.
order=1
fixediv=None
date0_inc=startdate-maxback
ok=False

onsprev,onsinc=getonsprevinc()

now=apiday()
while 1:
  data=getcases_raw(now,location="England")
  if 'Bad' not in data: break
  print("Can't get api data for %s. Backtracking to most recent usable date."%now)
  now-=1
completionhistory=0
while 1:
  completionhistory+=7
  data0=getcases_raw(now-completionhistory,location="England")
  if 'Bad' not in data0: break
print("Using api data as of",now,"and comparing to",now-completionhistory)
print("Variance multiplier for incidence and case counts =",odp)
print("Order",order)
print()
data0={Date(d):c for (d,c) in data0.items()}
data={Date(d):c for (d,c) in data.items()}
last=max(data)
# Correct recent incomplete entries, based on previous week
for d in Daterange(last-6,last+1):
  data[d]=round(data[d]*data[d-completionhistory]/data0[d-completionhistory])

date0_cases=Date(min(data))
if date0_cases>startdate: raise RuntimeError("Date %s not found in case data"%startdate)
cases=list(data.values())

enddate=max(onsprev[-1][1],onsinc[-1][1],date0_cases+len(cases))

n=enddate-startdate

def diffmat(n,order):
  A=np.zeros([n,n])
  if order==1:
    for i in range(n-1):
      A[i,i]+=1
      A[i,i+1]+=-1
      A[i+1,i]+=-1
      A[i+1,i+1]+=1
  elif order==2:
    for i in range(n-2):
      A[i,i]+=1
      A[i,i+1]+=-2
      A[i,i+2]+=1
      A[i+1,i]+=-2
      A[i+1,i+1]+=4
      A[i+1,i+2]+=-2
      A[i+2,i]+=1
      A[i+2,i+1]+=-2
      A[i+2,i+2]+=1
  else:
    raise RuntimeError("Unrecognised order %d"%order)
  return A

iv=10
A=iv*diffmat(n,order)
b=np.zeros(n)
c=0
# minimising quadratic form x^t.A.x-2b^t.x+c
# prob is proportional to exp(-(1/2)QF(x))

for date0,date1,targ,low,high in onsinc:
  targ/=scale
  var=((high-low)/scale/(2*1.96))**2*odp
  i0=date0-startdate
  i1=date1-startdate
  if i0>=0:
    A[i0:i1,i0:i1]+=1/((i1-i0)**2*var)
    b[i0:i1]+=targ/((i1-i0)*var)
    c+=targ**2/var

# Positivity kernel. Set to something simple pro tem
poskern=np.array([1]*9)
posback=7
k=len(poskern)
for date0,date1,targ,low,high in onsprev:
  targ/=scale
  var=((high-low)/scale/(2*1.96))**2*odp
  i0=date0-posback-startdate
  i1=date1-posback-startdate
  j0=date0-startdate
  j1=date1-startdate
  if i0>=0 and i1+k-1<=n:
    print(i0,i1,targ,sqrt(var))
    h=i1-i0
    poskern2=np.zeros(k+h-1)
    # We're assuming ONS positivity for a date range refers to the probabilty that a random person _on a random date within that range_ would test positive
    # (As opposed to, for example, the probability that a random person would test positive on some date within that range.)
    for i in range(h): poskern2[i:i+k]+=poskern/h
    A[i0:i1+k-1,i0:i1+k-1]+=np.outer(poskern2,poskern2)/var
    b[i0:i1+k-1]+=poskern2*targ/var
    c+=targ**2/var

    
xx=np.linalg.solve(A,b)
with open('temp%d'%order,'w') as fp:
  for (i,x) in enumerate(xx):
    print(startdate+i,x,file=fp)

sa

def getbackparams(back):
  d0_i=d0_c-back
  incdata=np.array(incidence[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
  vardata=np.array(varinc[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
  assert len(incdata)==n and len(vardata)==n
  var=((casedata/incdata)**2*vardata+casedata)*odp
  al=incdata*incdata/var
  b=incdata*casedata/var
  c=(casedata*casedata/var).sum()
  
# var[i] = (x[i]**2*V[incdata[i]]+V[casedata[i]])*odp
#       ~= ((case[i]/inc[i])**2*V[incdata[i]]+casedata[i])*odp
#
# First order version (CARs are similar over time):
#
# Q(x[], al[], iv) = sum_i (x[i]*incdata[i]-casedata[i])^2/var[i] + iv*sum_i (x[i+1]-x[i])^2
#                  = sum_i al[i].(x[i]-t[i])^2 + iv*sum_i (x[i+1]-x[i])^2, where t[i]=case[i]/inc[i], al[i]=inc[i]^2/var[i]
#                  = sum_{i,j} A_{i,j}x[i]x[j] - 2.sum_i b[i]x[i] + c
# Find iv by maximising \int_{x[]} exp(-(1/2)Q(x[], al[], iv)) / lim_{eps[]->0} (sum(eps[]))^{1/2}.\int_{x[]} exp(-(1/2)Q(x[], eps[], iv))
#                     = \int_{x[]} exp(-(1/2)Q(x[], al[], iv)) / iv^{-(1/2)(n-1)}
#
#
# Second order version (rate of change of CARs are similar over time):
#
# Q(x[], al[], iv) = sum_i (x[i]*incdata[i]-casedata[i])^2/var[i] + iv*sum_i (x[i+2]-2*x[i+1]+x[i])^2
#                  = sum_i al[i].(x[i]-t[i])^2 + iv*sum_i (x[i+2]-2*x[i+1]+x[i])^2, where t[i]=case[i]/inc[i], al[i]=inc[i]^2/var[i]
# Find iv by maximising \int_{x[]} exp(-(1/2)Q(x[], al[], iv)) / lim_{eps[]->0} (thing2(eps[])^{1/2}.\int_{x[]} exp(-(1/2)Q(x[], eps[], iv))
#                     = \int_{x[]} exp(-(1/2)Q(x[], al[], iv)) / iv^{-(1/2)(n-2)}
# 
# 
# (See gaussianrenormalisationcheck.py)
# Then use this iv to find x[] by minimising Q(x[], al[], iv)
def MLE_and_integrate(iv,al,b,c,order):
  n=len(b)
  A=iv*diffmat(n,order)
  for i in range(n): A[i,i]+=al[i]
  xx=np.linalg.solve(A,b)
  LL=(1/2)*(b@xx)-(1/2)*np.linalg.slogdet(A)[1]-(1/2)*c+(1/2)*(n-order)*log(iv)
  return xx,LL

def getcoeffs(back):
  d0_i=d0_c-back
  n=len(casedata)
  incdata=np.array(incidence[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
  vardata=np.array(varinc[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
  assert len(incdata)==n and len(vardata)==n
  tt=casedata/incdata
  ivar=1/((tt**2*vardata+casedata)*odp)
  for date in ignore:
    if (date-d0_c)%7==0:
      i=(date-d0_c)//7
      if i>=0 and i<n: ivar[i]=0
  al=incdata*incdata*ivar
  b=incdata*casedata*ivar
  c=(casedata*casedata*ivar).sum()
  return al,b,c,tt

# To make comparison fair, ensure same case days are predicted, regardless of back
# For all minback<=back<=maxback, ncasedays-back<=len(incidence)-(startdate-date0_inc)
ncasedays=len(incidence)-(startdate-date0_inc)+minback
date1_cases=date0_inc+len(incidence)+minback
totresid0=totresid1=totLL=0
for day in range(7):
  d0_c=startdate+(monday+day-startdate)%7
  casedata=np.array(cases[d0_c-date0_cases:date1_cases-date0_cases:7])
  n=len(casedata);assert n>=2
  best=(-1e30,)
  for back in range(minback,maxback+1):
    al,b,c,tt=getcoeffs(back)
    def NLL(logiv):
      xx,LL=MLE_and_integrate(exp(logiv),al,b,c,order)
      return -LL
    if fixediv==None:
      bounds=(-10,15)
      res=minimize(NLL,[0],bounds=[bounds],method="SLSQP")
      if not res.success: raise RuntimeError(res.message)
      if res.x<bounds[0]+1e-6 or res.x>bounds[1]-1e-6: print("Warning:",res.x,"touching bound")
      iv=exp(res.x)
      LL=-res.fun
    else:
      iv=fixediv
      xx,LL=MLE_and_integrate(iv,al,b,c,order)
    if LL>best[0]: best=[LL,back,iv]
  back=best[1]
  iv=best[2]
  al,b,c,tt=getcoeffs(back)
  xx,LL=MLE_and_integrate(iv,al,b,c,order)
  resid=sqrt((al*(xx-tt)**2).sum()/n)
  print("Day %d    back %d    iv %7.1f    average_CAR %6.3f    average_resid(sd units) %5.3f"%(day,back,iv,xx.sum()/n,resid))
  with open('temp.%d.%d.%d'%(order,day,odp),'w') as fp:
    for i in range(n): print(d0_c+i*7,xx[i],tt[i],sqrt(al[i])*(xx[i]-tt[i]),al[i],file=fp)
  totresid0+=n
  totresid1+=(al*(xx-tt)**2).sum()
  totLL+=LL

print("Overall average residual (in sd units)",sqrt(totresid1/totresid0))
print("Total log likelihood = %.1f"%totLL)
