import requests,sys
import numpy as np
from math import log,exp,sqrt
from stuff import *
from scipy.optimize import minimize
from random import normalvariate as normal,random,seed,randrange

# Comparing epiforecast incidence with case rates

startdate=Date("2021-01-01")
if len(sys.argv)>1: startdate=Date(sys.argv[1])
print("Start date =",startdate)

# population=55.98e6
monday=Date("2022-01-03")

minback=1
maxback=10
odp=7
# Order=1 if you think the prior for xx[] (CAR) is Brownian motion-like (in particular, Markov)
# Order=2 if you think the prior for xx[] is more like the integral of Brownian motion.
# I think xx[] does actually have "momentum" sometimes - there is a downwards (or upwards) trend that is predictive, so perhaps order=2 is better?
# Though it can dart around sometimes for random reasons, like Christmas.
order=2
fixediv=None
date0_inc=startdate-maxback
ok=False
incidence=[];varinc=[]
with open("Englandincidence_fromepiforecasts") as fp:
  for x in fp:
    y=x.strip().split()
    if y[0]==date0_inc: ok=True
    if y[0]>=date0_inc: incidence.append(float(y[1]));varinc.append(float(y[2])**2)
if not ok: raise RuntimeError("Date %s not found in incidence data"%date0_inc)

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

# Make synthetic data
if 0:
  s=randrange(1000000000)
  #s=714180712
  #s=161012188
  s=798489798
  print("Seed",s)
  seed(s)
  date0=min(date0_inc,date0_cases)
  date1=max(date0_inc+len(incidence),date0_cases+len(cases))
  N=date1-date0
  infections=np.zeros(N)
  lmean=10
  v=lmean
  for i in range(N):
    infections[i]=exp(v)
    v+=-0.001*(v-lmean)+0.1*normal(0,1)
  for i in range(len(incidence)):
    incidence[i]=infections[date0_inc+i-date0]
    varinc[i]=incidence[i]/odp
  mincar=0.025
  maxcar=0.9
  for day in range(7):
    day0=(monday-date0_cases+day)%7
    prev=np.zeros(order+1)
    prev[-1]=0.5# Initial CAR
    sd=0.05
    with open("tempsyncar.%d.%d.%d"%(order,day,odp),'w') as fp:
      for i in range(day0,len(cases),7):
        prev[0]=sd*normal(0,1)
        for j in range(order): prev[j+1]+=prev[j]
        car=min(max(prev[-1],mincar),maxcar)
        ni=infections[date0_cases+i-date0]
        c=car*(ni+sqrt(ni/odp)*normal(0,1))
        cases[i]=max(c+sqrt(c/odp)*normal(0,1),1)
        print(date0_cases+i,car,cases[i],file=fp)

with open("incidence_vs_adjcases","w") as fp:
  adjcases=weekdayadj(cases)
  adjcases=adjcases[startdate-date0_cases:]
  for i in range(max(len(incidence),len(adjcases))):
    print(startdate+i,end="",file=fp)
    if i<len(incidence): print("  %7d"%incidence[i],end="",file=fp)
    else: print("  %7s"%"-",end="",file=fp)
    if i<len(adjcases): print("  %7d"%adjcases[i],end="",file=fp)
    else: print("  %7s"%"-",end="",file=fp)
    print(file=fp)

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
  incdata=np.array(incidence[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
  vardata=np.array(varinc[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
  assert len(incdata)==n and len(vardata)==n
  var=((casedata/incdata)**2*vardata+casedata)*odp
  al=incdata*incdata/var
  b=incdata*casedata/var
  c=(casedata*casedata/var).sum()
  tt=casedata/incdata
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
