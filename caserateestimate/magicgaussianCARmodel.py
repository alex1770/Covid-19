import requests,sys
import numpy as np
from math import log,exp
from stuff import *
from scipy.optimize import minimize

# Comparing epiforecast incidence with case rates

startdate=Date("2021-01-01")
if len(sys.argv)>1: startdate=Date(sys.argv[1])
print("Start date =",startdate)

# population=55.98e6
monday=Date("2022-01-03")

minback=1
maxback=10
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
data0={Date(d):c for (d,c) in data0.items()}
data={Date(d):c for (d,c) in data.items()}
last=max(data)

# Correct recent incomplete entries, based on previous week
for d in Daterange(last-6,last+1):
  data[d]=round(data[d]*data[d-completionhistory]/data0[d-completionhistory])

date0_cases=Date(min(data))
if date0_cases>startdate: raise RuntimeError("Date %s not found in case data"%startdate)
cases=list(data.values())

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

# To make comparison fair, ensure same case days are predicted, regardless of back
# For all minback<=back<=maxback, ncasedays-back<=len(incidence)-(startdate-date0_inc)
ncasedays=len(incidence)-(startdate-date0_inc)+minback
date1_cases=date0_inc+len(incidence)+minback
for day in range(7):
  d0_c=startdate+(monday+day-startdate)%7
  casedata=np.array(cases[d0_c-date0_cases:date1_cases-date0_cases:7])
  n=len(casedata);assert n>=2
  best=(1e30,)
  for back in range(minback,maxback+1):
    d0_i=d0_c-back
    incdata=np.array(incidence[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
    vardata=np.array(varinc[d0_i-date0_inc:d0_i-date0_inc+7*n:7])
    assert len(incdata)==n and len(vardata)==n
    var=(casedata/incdata)**2*vardata+casedata
    al=incdata*incdata/var
    b=incdata*casedata/var
    c=(casedata*casedata/var).sum()
    def NLL(logiv):
      # Perhaps add overdispersion parameter
      # var[i] = x[i]**2*V[incdata[i]]+V[casedata[i]]
      #       ~= (case[i]/inc[i])**2*V[incdata[i]]+casedata[i]
      # Q(x[], al[], iv) = sum_i { (x[i]*incdata[i]-casedata[i])^2/var[i] + iv*sum_i (x[i+1]-x[i])^2 }
      #                  = sum_i { inc[i]^2/var[i]*(x[i]-case[i]/inc[i])^2 + iv*sum_i (x[i+1]-x[i])^2 }
      #                  = sum_i { al[i].(x[i]-t[i])^2 + iv*sum_i (x[i+1]-x[i])^2 }, where t[i]=case[i]/inc[i], al[i]=inc[i]^2/var[i]
      # Find iv by maximising \int_{x[]} exp(-(1/2)Q(x[], al[], iv)) / lim_{eps[]->0} (sum(eps[]))^{1/2}.\int_{x[]} exp(-(1/2)Q(x[], eps[], iv))
      #                     = \int_{x[]} exp(-(1/2)Q(x[], al[], iv)) / iv^{-(1/2)(n-1)}
      # (See gaussianrenormalisationcheck.py)
      # Then use this iv to find x[] by minimising Q(x[], al[], iv)
      iv=exp(logiv)
      A=np.zeros([n,n])
      for i in range(n):
        A[i,i]+=al[i]
      for i in range(n-1):
        A[i,i]+=iv
        A[i+1,i+1]+=iv
        A[i,i+1]-=iv
        A[i+1,i]-=iv
      xx=np.linalg.solve(A,b)
      LL=(1/2)*(b@xx)-(1/2)*np.linalg.slogdet(A)[1]-(1/2)*c+(1/2)*(n-1)*log(iv)
      return -LL
    
    res=minimize(NLL,[0],bounds=[(-10,15)],method="SLSQP")#,options=minopts)
    if not res.success: raise RuntimeError(res.message)
    iv=exp(res.x)
    print("%d %2d %8.3f %9.3f"%(day,back,iv,res.fun))
    if res.fun<best[0]: best=[res.fun,back]
  print(day,best[1])
  print()
  
