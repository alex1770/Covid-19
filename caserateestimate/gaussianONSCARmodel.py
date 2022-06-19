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

ignore={Date(date) for date in ignore}

# population=55.98e6
monday=Date("2022-01-03")
scale=1e5# Units of this number of people to keep things nicely conditioned

np.set_printoptions(precision=3,suppress=True,linewidth=200)

minback=1
maxback=10
# odp_ons is an overdispersion parameter. Corrects ONS variance estimates which I think are significantly too low.
odp_ons=4         # Inverse coupling of incidence to ONS prevalence
odp_casedata=1    # Inverse coupling of incidence to case data
iv_inc=10         # Coupling of incidence to iteself
iv_car=30         # Coupling of inverse-CAR to itself
# Order=1 if you think the prior is exp(Brownian motion)-like (in particular, Markov)
# Order=2 if you think the prior is more like exp(integral of Brownian motion).
order=2

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
print("Variance multiplier for incidence and case counts =",odp_ons)
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

def difflogmat(n,order,xx0):
  invxx0=1/xx0
  logxx0=np.log(xx0)
  A=np.zeros([n,n])
  b=np.zeros(n)
  c=0
  v=[1]
  for i in range(order): v=np.pad(v,(0,1),mode="constant")-np.pad(v,(1,0),mode="constant")
  for i in range(n-order):
    w=v*invxx0[i:i+order+1]
    l=v@logxx0[i:i+order+1]
    A[i:i+order+1,i:i+order+1]+=np.outer(w,w)
    b[i:i+order+1]+=-w*l
    c+=l*l
  return A,b,c


# Positivity kernel. https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-021-01982-x, Fig 3b
# P(PCR detection at x days after infection)
# Was going to adjust this by a self-determined slowly-moving constant factor, but after looking at it,
# I believe there isn't good evidence for the constant being different from 1.

poskern=np.array([x*exp(-x/5)*0.45 for x in range(19,0,-1)])
numk=len(poskern)

casedata=np.array(cases[startdate-date0_cases:])/scale
ncases=len(casedata);assert ncases<=N
back=[5]*7# Pro tem

# 2N variables:
#    i       I[i] = incidence variable 0<=i<N
#  N+i       c[i] = CAR variable       0<=i<N

def savevars(xx,name="temp"):
  with open(name+"_inc",'w') as fp:
    for i in range(N):
      print(startdate+i,xx[i],file=fp)
  with open(name+"_CARcases",'w') as fp:
    for j in range(N):
      if startdate+j not in ignore:
        day=(startdate+j-monday)%7
        i=j-back[day]
        if i>=0:
          print(startdate+i,1/xx[N+j],xx[N+j]*casedata[j],file=fp)
        
  

def calibrateprevalencetoincidence():
  # Use ONS incidence as a one-off to calibrate use of ONS prevalence, in particular
  # poskern. The ratio of (ONSincidence convolved with poskernel) to ONSprevalence varies
  # somewhat around 1, though not in an obvious pattern that depends on variants. Therefore
  # it seems most likely (or at least, no evidence against the idea) that this variation is
  # part of ONS incidence being messed up, which is a known thing, which mean we discard ONS
  # incidence for the next stage, and anchor to (calibrated) ONS prevalence.

  iv_inc=100
  A=iv_inc*diffmat(N,order)
  b=np.zeros(N)
  c=0
  # minimising quadratic form x^t.A.x-2b^t.x+c
  # prob is proportional to exp(-(1/2)QF(x))
  adjinc=0
  for date0,date1,targ,low,high in onsinc:
    targ/=scale
    var=((high-low)/scale/(2*1.96))**2*odp_ons
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
      var=((high-low)/scale/(2*1.96))**2*odp_ons
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

def initialguess():
  
  A=iv_inc*diffmat(N,order)
  b=np.zeros(N)
  c=0
  
  for date0,date1,targ,low,high in onsprev:
    targ/=scale
    var=((high-low)/scale/(2*1.96))**2*odp_ons
    i0=date0-numk-startdate
    i1=date1-numk-startdate
    if i0>=0 and i1+numk-1<=N:
      h=i1-i0
      poskern2=np.zeros(numk+h-1)
      # We're assuming ONS positivity for a date range refers to the probabilty that a random person _on a random date within that range_ would test positive
      # (As opposed to, for example, the probability that a random person would test positive on some date within that range.)
      for i in range(h): poskern2[i:i+numk]+=poskern/h
      A[i0:i1+numk-1,i0:i1+numk-1]+=np.outer(poskern2,poskern2)/var
      b[i0:i1+numk-1]+=poskern2*targ/var
      c+=targ**2/var
  
  # Initial incidence guess
  inc0=np.maximum(np.linalg.solve(A,b),0.01)
  
  A=iv_car*diffmat(N,order)
  b=np.zeros(N)
  c=0

  # Link cases to incidence using al[i].(I[i]-c[i'].casedata[i'])^2
  #                               al[i]=(c[i']^2.V[casedata[i']])^{-1}
  #                                   ~=((I[i]/casedata[i'])^2.casedata[i'].odp_casedata)^{-1}
  #                                    =(odp_casedata.I[i']^2/casedata[i'])^{-1}
  for j in range(ncases):
    if startdate+j not in ignore:
      day=(startdate+j-monday)%7
      i=j-back[day]
      if i>=0:
        al=casedata[j]/(odp_casedata*inc0[i]**2)
        A[j,j]+=al*casedata[j]**2
        b[j]+=al*casedata[j]*inc0[i]
        #print(startdate+j,inc0[i]/casedata[j])
  
  # Initial inverse-CAR guess
  c0=np.linalg.solve(A,b)
  
  xx=np.zeros(N*2)
  xx[:N]=inc0
  xx[N:]=c0
  return xx

xx=initialguess()
savevars(xx,"tempinit")

for it in range(100):
  print("Iteration",it)
  A=np.zeros([N*2,N*2])
  b=np.zeros(N*2)
  c=0

  A_i,b_i,c_i=difflogmat(N,order,xx[:N])
  A[:N,:N]+=iv_inc*A_i
  b[:N]+=iv_inc*b_i
  c+=iv_inc*c_i
  
  # Link incidence to ONSprevalence
  for date0,date1,targ,low,high in onsprev:
    targ/=scale
    var=((high-low)/scale/(2*1.96))**2*odp_ons
    i0=date0-numk-startdate
    i1=date1-numk-startdate
    if i0>=0 and i1+numk-1<=N:
      h=i1-i0
      poskern2=np.zeros(numk+h-1)
      # We're assuming ONS positivity for a date range refers to the probabilty that a random person _on a random date within that range_ would test positive
      # (As opposed to, for example, the probability that a random person would test positive on some date within that range.)
      for i in range(h): poskern2[i:i+numk]+=poskern/h
      A[i0:i1+numk-1,i0:i1+numk-1]+=np.outer(poskern2,poskern2)/var
      b[i0:i1+numk-1]+=poskern2*targ/var
      c+=targ**2/var
  
  
  
  # Add in diff constraints for CAR variables which relate same days of week to each other
  # (d isn't necessarily equal to the day of the week)
  for d in range(7):
    n=(enddate-(startdate+d)+6)//7
    A_c,b_c,c_c=difflogmat(n,order,xx[N+d::7])
    c+=iv_car*c_c
    for i in range(n):
      b[N+d+i*7]+=iv_car*b_c[i]
      for j in range(n):
        A[N+d+i*7,N+d+j*7]+=iv_car*A_c[i,j]
  
  # Terms al[i].(I[i]-c[i'].casedata[i'])^2 correspond to CAR error at incidence i (casedata i')
  #       al[i]=(V[I[i]-c[i'].casedata[i']])^{-1}
  #           ~=(odp_casedata.I0[i]^2/casedata[i'])^{-1}
  #             where I0[i] is some approximation to the incidence
  al=np.zeros(N)
  for j in range(ncases):
    if startdate+j not in ignore:
      day=(startdate+j-monday)%7
      i=j-back[day]
      if i>=0:
        al[i]=casedata[j]/(odp_casedata*xx[i]**2)
        A[i,i]+=al[i]
        A[i,N+j]-=al[i]*casedata[j]
        A[N+j,i]-=al[i]*casedata[j]
        A[N+j,N+j]+=al[i]*casedata[j]**2

  xx0=xx
  xx=np.linalg.solve(A,b)
  if np.abs(xx/xx0-1).max()<1e-6: break

savevars(xx)
