import requests,sys
import numpy as np
from math import log,exp,sqrt
from stuff import *
from scipy.optimize import minimize
from random import normalvariate as normal,random,seed,randrange
from parseonsprevinc import getonsprevinc

startdate=Date("2020-06-01")
if len(sys.argv)>1: startdate=Date(sys.argv[1])
print("Start date =",startdate)

# Dates when case data is unreliable and is to be ignored (e.g., Christmas)
ignore={"2020-06-11", "2020-06-28", "2020-07-04", "2020-07-12", "2020-07-19", "2020-07-31", "2020-08-02", "2020-08-09", "2020-08-10", "2020-08-11", "2020-08-15",
        "2020-08-16", "2020-08-17", "2020-08-31", "2020-09-01", "2020-09-13", "2020-09-21", "2020-11-02", "2020-12-19", "2020-12-25", "2020-12-26",
        "2020-12-27", "2020-12-28", "2020-12-31", "2021-01-01", "2021-01-23", "2021-01-24", "2021-03-21", "2021-03-28", "2021-04-02", "2021-04-05",
        "2021-04-18", "2021-05-09", "2021-05-31", "2021-08-30", "2021-08-31", "2021-12-25", "2021-12-29", "2022-01-01", "2022-01-04", "2022-02-18",
        "2022-12-25", "2022-12-26", "2022-12-27"}
ignore={Date(date) for date in ignore}

# population=55.98e6
monday=Date("2022-01-03")
scale=1e5# Units of this number of people to keep things nicely conditioned

np.set_printoptions(precision=3,suppress=True,linewidth=200)

minback=1
maxback=10
inc_ons=1/4    # Coupling of incidence to ONS prevalence (less than 1 means we think ONS confidence intervals are too narrow)
inc_case=1     # Coupling of incidence to case data (less than 1 means we think case data is "overdispersed" with a variance bigger than the count)
inc_inc=10000    # Coupling of incidence to iteself
car_car=3000     # Coupling of inverse-CAR to itself
# Order=1 if you think the prior is exp(Brownian motion)-like (in particular, Markov)
# Order=2 if you think the prior is more like exp(integral of Brownian motion).
order=2

def meanvar(sample):
  mu=sum(sample)/len(sample)
  var=sum((x-mu)**2 for x in sample)/(len(sample)-1)
  return mu,var

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
#poskern=np.array([10,0,0,0,0,0,0,0,0,0])
#poskern=np.zeros(21)+0.53
#poskern=np.zeros(50)+0.2
numk=len(poskern)

# 2N variables:
#    i       I[i] = incidence variable 0<=i<N
#  N+i       c[i] = CAR variable       0<=i<N

def savevars(N,casedata,back,xx,name="temp"):
  with open(name+"_incidence",'w') as fp:
    for i in range(N):
      print(startdate+i,"%9.1f"%(xx[i]*scale),file=fp)
  with open(name+"_CARcases",'w') as fp:
    for j in range(N):
      if startdate+j not in ignore:
        day=(startdate+j-monday)%7
        i=j-back[day]
        if i>=0:
          print(startdate+i,"%7.4f %9.1f"%(1/xx[N+j],xx[N+j]*casedata[j]*scale),file=fp)



def calibrateprevalencetoincidence():
  # Use ONS incidence as a one-off to calibrate use of ONS prevalence, in particular
  # poskern. The ratio of (ONSincidence convolved with poskernel) to ONSprevalence varies
  # somewhat around 1, though not in an obvious pattern that depends on variants. Therefore
  # it seems most likely (or at least, no evidence against the idea) that this variation is
  # part of ONS incidence being messed up, which is a known thing, which mean we discard ONS
  # incidence for the next stage, and anchor to (calibrated) ONS prevalence.

  inc_inc=100
  A=inc_inc*diffmat(N,order)
  b=np.zeros(N)
  c=0
  # minimising quadratic form x^t.A.x-2b^t.x+c
  # prob is proportional to exp(-(1/2)QF(x))
  adjinc=0
  for date0,date1,targ,low,high in onsinc:
    targ/=scale
    var=((high-low)/scale/(2*1.96))**2/inc_ons
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
      var=((high-low)/scale/(2*1.96))**2/inc_ons
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

def getcaseoutliers():
  first=1
  scd=np.log(xx[N:]*casedata)
  for j in range(1,N-1):
    if startdate+j not in ignore:
      l=[]
      for k in [j-1,j+1]:
        if startdate+k not in ignore: l.append(scd[k])
      if len(l)>0:
        o=scd[j]-sum(l)/len(l)
        if abs(o)>0.2:
          if first: print("Case outliers not already ignored (sample days shown):");first=0
          print("Case outlier:",startdate+j,"%6.3f"%o)

def initialguess(N,onsprev,casedata,back):
  
  A=inc_inc*diffmat(N,order)
  b=np.zeros(N)
  c=0
  
  for date0,date1,targ,low,high in onsprev:
    targ/=scale
    var=((high-low)/scale/(2*1.96))**2/inc_ons
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
  
  A=car_car*diffmat(N,order)
  b=np.zeros(N)
  c=0

  # Link cases to incidence using al[i].(I[i]-c[i'].casedata[i'])^2
  #                               al[i]=(c[i']^2.V[casedata[i']])^{-1}
  #                                   ~=((I[i]/casedata[i'])^2.(casedata[i'].inc_case))^{-1}
  #                                    =(I[i']^2/(casedata[i'].inc_case))^{-1}
  for j in range(len(casedata)):
    if startdate+j not in ignore:
      day=(startdate+j-monday)%7
      i=j-back[day]
      if i>=0:
        al=casedata[j]*inc_case/inc0[i]**2
        A[j,j]+=al*casedata[j]**2
        b[j]+=al*casedata[j]*inc0[i]
        #print(startdate+j,inc0[i]/casedata[j])
  
  # Initial inverse-CAR guess
  c0=np.maximum(np.linalg.solve(A,b),0.01)
  
  xx=np.zeros(N*2)
  xx[:N]=inc0
  xx[N:]=c0
  return xx

def getest(now=apiday(),prlev=0):
  now=Date(now)
  onsprev,onsinc=getonsprevinc(maxdate=now)
  while 1:
    data=getcases_raw(now,location="England")
    if 'Bad' not in data: break
    if prlev>=2: print("Can't get api data for %s. Backtracking to most recent usable date."%now)
    now-=1
  ncompsamples=6
  completionlist=[]
  completionhistory=0
  while len(completionlist)<ncompsamples:
    completionhistory+=7
    data0=getcases_raw(now-completionhistory,location="England")
    if 'Bad' not in data0: completionlist.append((completionhistory,data0))
  if prlev>=2:
    print("Using api data as of",now,"and comparing to",' '.join(str(now-ch) for (ch,data0) in completionlist))
    print("Order",order)
    print()
  last=max(data)
  # Correct recent incomplete entries, based on the same day in previous weeks
  for d in Daterange(last-6,last+1):
    sample=[log(data[d-ch]/data0[d-ch]) for (ch,data0) in completionlist]
    mu,var=meanvar(sample)
    #print(d,"%7.3f %7.3f"%(mu,sqrt(var)))
    data[d]=round(data[d]*exp(mu))
  
  date0_cases=Date(min(data))
  if date0_cases>startdate: raise RuntimeError("Date %s not found in case data"%startdate)
  cases=list(data.values())
  #for i in range(14): cases[i-14]*=exp(-i*0.05)

  enddate=max(onsprev[-1][1],onsinc[-1][1],date0_cases+len(cases))
  # Possibly extend this to extrapolate
  N=enddate-startdate
  
  casedata=np.array(cases[startdate-date0_cases:])/scale
  assert len(casedata)<=N
  back=[5]*7# Pro tem
  
  
  xx=initialguess(N,onsprev,casedata,back)
  savevars(N,casedata,back,xx,"tempinit")
  
  for it in range(20):
    if prlev>=2: print("Iteration",it)
    A=np.zeros([N*2,N*2])
    b=np.zeros(N*2)
    c=0
  
    A_i,b_i,c_i=difflogmat(N,order,xx[:N])
    A[:N,:N]+=inc_inc*A_i
    b[:N]+=inc_inc*b_i
    c+=inc_inc*c_i
    
    # Link incidence to ONSprevalence
    for date0,date1,targ,low,high in onsprev:
      targ/=scale
      var=((high-low)/scale/(2*1.96))**2/inc_ons
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
      n=(N-d+6)//7
      A_c,b_c,c_c=difflogmat(n,order,xx[N+d::7])
      c+=car_car*c_c
      for i in range(n):
        b[N+d+i*7]+=car_car*b_c[i]
        for j in range(n):
          A[N+d+i*7,N+d+j*7]+=car_car*A_c[i,j]
    
    # Terms al[i].(I[i]-c[i'].casedata[i'])^2 correspond to CAR error at incidence i (casedata i')
    #       al[i]=(V[I[i]-c[i'].casedata[i']])^{-1}
    #           ~=(I0[i]^2/(casedata[i'].inc_case))^{-1}
    #             where I0[i] is some approximation to the incidence
    al=np.zeros(N)
    for j in range(len(casedata)):
      if startdate+j not in ignore:
        day=(startdate+j-monday)%7
        i=j-back[day]
        if i>=0:
          al[i]=casedata[j]*inc_case/xx[i]**2
          A[i,i]+=al[i]
          A[i,N+j]-=al[i]*casedata[j]
          A[N+j,i]-=al[i]*casedata[j]
          A[N+j,N+j]+=al[i]*casedata[j]**2
  
    xx0=xx
    xx=np.maximum(np.linalg.solve(A,b),0.01)
    if np.abs(xx/xx0-1).max()<1e-3: break

  savevars(N,casedata,back,xx,"England")
  return casedata,xx
  
  
  
casedata0,xx0=getest()
#getcaseoutliers()
N=xx0.shape[0]//2


while 1:
  t=random()*2-1
  u=random()*2-1
  v=random()*2-1
  w=random()*2-1
  inc_inc=exp(t)*135
  car_car=exp(u)*41
  inc_ons=exp(v)*0.2
  inc_case=exp(w)*1.0
    
  now=apiday()
  numcheck=30
  chrange=7
  err=0
  for ch in range(numcheck):
    casedata,xx=getest(now-(ch+1)*chrange,prlev=0)
    i0,i1=N-(ch+2)*chrange,N-(ch+1)*chrange
    #print(xx0[i0:i1])
    #print(xx[i0:i1])
    #print()
    err+=(np.log(xx[i0:i1]/xx0[i0:i1])**2).sum()
    xx=xx0
  
  print("%12g %12g %12g %12g     %7.3f"%(inc_ons,inc_case,inc_inc,car_car,sqrt(err/(numcheck*chrange))))
  sys.stdout.flush()
        
