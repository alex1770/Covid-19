import requests,sys
import numpy as np
from math import log,exp,sqrt
from stuff import *
from scipy.optimize import minimize
from random import normalvariate as normal,random,seed,randrange
from parseonsprevinc import getonsprevinc

startdate=Date("2021-05-31")
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
odp=4# Overdispersion parameter. Corrects ONS variance estimates which I think are significantly too low.
odp_casedata=1
# Order=1 if you think the prior for xx[] is Brownian motion-like (in particular, Markov)
# Order=2 if you think the prior for xx[] is more like the integral of Brownian motion.
order=2
fixediv=None
date0_inc=startdate-maxback
iv_inc=100
iv_car=0.01

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

# Possibly extend this to extrapolate
N=enddate-startdate

# Returns matrix that calculates 'order' order differences
# E.g., for order=1, QF(x)=sum_i (x_i-x_{i+1})^2
# [  1 -1  0  0 ... ]
# [ -1  2 -1  0 ... ]
# [  0 -1  2 -1 ... ]
# [  0  0 -1  2 ... ]
# ...
def diffmat(n,order):
  A=np.zeros([n,n])
  v=[1]
  for i in range(order): v=np.pad(v,(0,1),mode="constant")-np.pad(v,(1,0),mode="constant")
  vv=np.outer(v,v)
  for i in range(n-order): A[i:i+order+1,i:i+order+1]+=vv
  return A


# Positivity kernel. https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-021-01982-x, Fig 3b
# P(PCR detection at x days after infection)
# Was going to adjust this by a self-determined slowly-moving constant factor, but after looking at it,
# I believe there isn't good evidence for the constant being different from 1.

poskern=np.array([x*exp(-x/5)*0.45 for x in range(19,0,-1)])
numk=len(poskern)

# Use ONS incidence as a one-off to calibrate use of ONS prevalence, in particular
# poskern. The ratio of (ONSincidence convolved with poskernel) to ONSprevalence varies
# somewhat around 1, though not in an obvious pattern that depends on variants. Therefore
# it seems most likely (or at least, no evidence against the idea) that this variation is
# part of ONS incidence being messed up, which is a known thing, which mean we discard ONS
# incidence for the next stage, and anchor to (calibrated) ONS prevalence.

if 0:
  iv_inc=100
  A=iv_inc*diffmat(N,order)
  b=np.zeros(N)
  c=0
  # minimising quadratic form x^t.A.x-2b^t.x+c
  # prob is proportional to exp(-(1/2)QF(x))
  adjinc=0
  for date0,date1,targ,low,high in onsinc:
    targ/=scale
    var=((high-low)/scale/(2*1.96))**2*odp
    i0=date0-adjinc-startdate
    i1=date1-adjinc-startdate
    if i0>=0 and i1<=N:
      A[i0:i1,i0:i1]+=1/((i1-i0)**2*var)
      b[i0:i1]+=targ/((i1-i0)*var)
      c+=targ**2/var
  
  inc_v0=np.linalg.solve(A,b)
  for adjinc in range(-10,21):
    s0=s1=s2=0
    for date0,date1,targ,low,high in onsprev:
      targ/=scale
      var=((high-low)/scale/(2*1.96))**2*odp
      i0=date0-numk+adjinc-startdate
      i1=date1-numk+adjinc-startdate
      if i0>=0 and i1+numk-1<=N:
        h=i1-i0
        poskern2=np.zeros(numk+h-1)
        # We're assuming ONS positivity for a date range refers to the probabilty that a random person _on a random date within that range_ would test positive
        # (As opposed to, for example, the probability that a random person would test positive on some date within that range.)
        for i in range(h): poskern2[i:i+numk]+=poskern/h
        pprev=inc_v0[i0:i1+numk-1]@poskern2
        x=log(targ/pprev)
        s0+=1;s1+=x;s2+=x*x
        #print(date0,pprev,targ,x)
    sd=sqrt((s2-s1**2/s0)/(s0-1))
    print(adjinc,sd)
  qwe

A=iv_inc*diffmat(N,order)
b=np.zeros(N)
c=0

for date0,date1,targ,low,high in onsprev:
  targ/=scale
  var=((high-low)/scale/(2*1.96))**2*odp
  i0=date0-numk-startdate
  i1=date1-numk-startdate
  if i0>=0 and i1+numk-1<=N:
    #print(i0,i1,targ,sqrt(var))
    h=i1-i0
    poskern2=np.zeros(numk+h-1)
    # We're assuming ONS positivity for a date range refers to the probabilty that a random person _on a random date within that range_ would test positive
    # (As opposed to, for example, the probability that a random person would test positive on some date within that range.)
    for i in range(h): poskern2[i:i+numk]+=poskern/h
    A[i0:i1+numk-1,i0:i1+numk-1]+=np.outer(poskern2,poskern2)/var
    b[i0:i1+numk-1]+=poskern2*targ/var
    c+=targ**2/var

inc0=np.linalg.solve(A,b)
with open("temp",'w') as fp:
  for (i,x) in enumerate(inc0):
    print(startdate+i,x,file=fp)

casedata=np.array(cases[startdate-date0_cases:])/scale
ncases=len(casedata);assert ncases<=N
# Pad to incorporate variables c(t) rerpesenting 1/CAR(t)
# In comments, denote i       I[i] = incidence variable 0<=i<N
#                     N+i  c[i] = CAR variable       0<=i<ncases
A=np.pad(A,(0,ncases),mode="constant")
b=np.pad(b,(0,ncases),mode="constant")

# Add in diff constraints for CAR variables which relate same days of week to each other
# (d isn't necessarily equal to the day of the week)
for d in range(7):
  n=(enddate-(startdate+d)+6)//7
  a=diffmat(n,order)
  for i in range(n):
    for j in range(n):
      A[N+d+i*7,N+d+j*7]+=iv_car*a[i,j]

# Terms al[i].(I[i]-c[i'].casedata[i'])^2 correspond to CAR error at incidence i (casedata i')
#       al[i]=V[I[i]-c[i'].casedata[i']]^{-1}
#           ~=(odp_casedata*I0[i]^2/casedata[i'])^{-1}
#             where I0[i] is some approximation to the incidence
back=[5]*7# Pro tem
al=np.zeros(N)
for j in range(ncases):
  if startdate+j not in ignore:
    day=(startdate+j-monday)%7
    i=j-back[day]
    if i>=0:
      al[i]=casedata[j]/(odp_casedata*inc0[i]**2)
      A[i,i]+=al[i]
      A[i,N+j]-=al[i]*casedata[j]
      A[N+j,i]-=al[i]*casedata[j]
      A[N+j,N+j]+=al[i]*casedata[j]**2
    

xx=np.linalg.solve(A,b)
with open("temp2",'w') as fp:
  for i in range(N):
    print(startdate+i,xx[i],file=fp)
with open("temp3",'w') as fp:
  for j in range(N):
    day=(startdate+j-monday)%7
    i=j-back[day]
    if i>=0:
      print(startdate+i,xx[N+j],xx[N+j]*casedata[j],file=fp)
asd


# Q(x[], al[], iv) = sum_i (x[i]*incdata[i]-casedata[i])^2/var[i] + iv*sum_i (x[i+1]-x[i])^2
#                  = sum_i al[i].(x[i]-t[i])^2 + iv*sum_i (x[i+1]-x[i])^2, where t[i]=case[i]/inc[i], al[i]=inc[i]^2/var[i]
#                  = sum_{i,j} A_{i,j}x[i]x[j] - 2.sum_i b[i]x[i] + c


d0_i=d0_c-back
incdata=np.array(incidence[d0_i-date0_inc:d0_i-date0_inc+7*N:7])
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
