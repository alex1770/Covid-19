import requests,sys
import numpy as np
from math import log,exp,sqrt
from stuff import *
from scipy.stats import multivariate_normal as mvn
from random import normalvariate as normal,random,seed,randrange
from parseonsprevinc import getonsprevinc,getdailyprevalence

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
back=[5]*7# Pro tem
inc_ons=exp(-3)     # Coupling of incidence to ONS prevalence (fixed) (ONS confidence intervals have now been adjusted to be sensible, but they are not independent)
inc_case=exp(3.0)   # Coupling of incidence to case data (less than 1 means we think case data is "overdispersed" with a variance bigger than the count)
inc_inc=exp(10.9)   # Coupling of incidence to iteself
car_car=exp(6.0)    # Coupling of inverse-CAR to itself week-by-week
car_car_d=exp(0.0)  # Coupling of inverse-CAR to itself day-by-day
# Order=1 if you think the prior is exp(Brownian motion)-like (in particular, Markov)
# Order=2 if you think the prior is more like exp(integral of Brownian motion) (has "momentum")
# etc
order=2

def rnd(): return random()*2-1

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

def savevars(N,casedata,back,xx,low=None,high=None,name="temp"):
  with open(name+"_incidence",'w') as fp:
    for i in range(N):
      print(startdate+i,"%9.1f"%(xx[i]*scale),end="",file=fp)
      if low is not None: print("  %9.1f"%(low[i]*scale),end="",file=fp)
      else: print("          -",end="",file=fp)
      if high is not None: print("  %9.1f"%(high[i]*scale),end="",file=fp)
      else: print("          -",end="",file=fp)
      print(file=fp)
  with open(name+"_CARcases",'w') as fp:
    for j in range(N):
      if startdate+j not in ignore:
        day=(startdate+j-monday)%7
        i=j-back[day]
        if i>=0:
          print(startdate+i,"%7.4f"%(1/xx[N+j]),end="",file=fp)
          if high is not None: print("  %7.4f"%(1/high[N+j]),end="",file=fp)
          else: print("        -",end="",file=fp)
          if low is not None: print("  %7.4f"%(1/low[N+j]),end="",file=fp)
          else: print("        -",end="",file=fp)
          if j<len(casedata): print("  %9.1f"%(xx[N+j]*casedata[j]*scale),end="",file=fp)
          else: print("          -",end="",file=fp)
          if j<len(casedata) and low is not None: print("  %7.4f"%(low[N+j]*casedata[j]*scale),end="",file=fp)
          else: print("        -",end="",file=fp)
          if j<len(casedata) and high is not None: print("  %7.4f"%(high[N+j]*casedata[j]*scale),end="",file=fp)
          else: print("        -",end="",file=fp)
          print(file=fp)



def calibrateprevalencetoincidence():
  # Not currently used, and not been updated since we went over to using daily ONS prevalence
  # 
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
        #print(date0,pprev,trg,x)
    sd=sqrt((s2-s1**2/s0)/(s0-1))
    print(adjinc,sd)

def getcaseoutliers(casedata,N):
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
  
  for date,targ,sd in onsprev:
    targ/=scale
    sd/=scale
    var=sd**2/inc_ons
    i=date-numk-startdate
    if i>=0 and i+numk<=N:
      A[i:i+numk,i:i+numk]+=np.outer(poskern,poskern)/var
      b[i:i+numk]+=poskern*targ/var
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
  c0=np.maximum(np.linalg.solve(A,b),1)
  
  xx=np.zeros(N*2)
  xx[:N]=inc0
  xx[N:]=c0
  return xx

def getextdata(enddate=apiday(),prlev=0):
  enddate=Date(enddate)
  onsprev=getdailyprevalence(maxdate=enddate)
  apireq=enddate
  while 1:
    data=getcases_raw(apireq,location="England")
    if 'Bad' not in data: break
    if prlev>=2: print("Can't get api data for %s. Backtracking to most recent usable date."%apireq)
    apireq-=1
  ncompsamples=6
  completionlist=[]
  completionhistory=0
  while len(completionlist)<ncompsamples:
    completionhistory+=7
    data0=getcases_raw(apireq-completionhistory,location="England")
    if 'Bad' not in data0: completionlist.append((completionhistory,data0))
  if prlev>=2:
    print("Using api data as of",apireq,"and comparing to",' '.join(str(apireq-ch) for (ch,data0) in completionlist))
    print("Order",order)
    print()
  last=max(data)
  # Correct recent incomplete entries, based on the same day in previous weeks
  for d in Daterange(last-6,last+1):
    sample=[log(data[d-ch]/data0[d-ch]) for (ch,data0) in completionlist]
    mu,var=meanvar(sample)
    # To do: use var
    #print(d,"%7.3f %7.3f"%(mu,sqrt(var)))
    data[d]=round(data[d]*exp(mu))
  
  date0_cases=Date(min(data))
  if date0_cases>startdate: raise RuntimeError("Date %s not found in case data"%startdate)
  cases=list(data.values())
  #for i in range(14): cases[i-14]*=exp(-i*0.05)
  #for i in range(14):
  #  pr=onsprev[i-14]
  #  onsprev[i-14]=(pr[0],pr[1]*exp(-i*0.05),pr[2])

  N=enddate-startdate
  
  casedata=np.array(cases[startdate-date0_cases:])/scale
  assert len(casedata)<=N
  
  return N,casedata,onsprev

# Returns two quadratic(etc) forms:
#   A0,b0,c0 <-> Prior distribution on internal variables x_i (incidence), a_i (inverse CAR)
#   A1,b1,c1 <-> Distribution of data (casecounts, ONS prevalence) given internal variables
def getqform(N,casedata,onsprev,xx0):
  
  A0=np.zeros([N*2,N*2])
  b0=np.zeros(N*2)
  c0=0
  A1=np.zeros([N*2,N*2])
  b1=np.zeros(N*2)
  c1=0
  numprev=numcases=0

  A_i,b_i,c_i=difflogmat(N,order,xx0[:N])
  A0[:N,:N]+=inc_inc*A_i
  b0[:N]+=inc_inc*b_i
  c0+=inc_inc*c_i
  
  # Link incidence to ONSprevalence
  for date,targ,sd in onsprev:
    targ/=scale
    sd/=scale
    var=sd**2/inc_ons
    i=date-numk-startdate
    if i>=0 and i+numk<=N:
      A1[i:i+numk,i:i+numk]+=np.outer(poskern,poskern)/var
      b1[i:i+numk]+=poskern*targ/var
      c1+=targ**2/var
      numprev+=1
  c1-=numprev*log(inc_ons)
  
  # Add in diff constraints for CAR variables which relate same days of week to each other
  # (d isn't necessarily equal to the day of the week)
  for d in range(7):
    n=(N-d+6)//7
    A_c,b_c,c_c=difflogmat(n,order,xx0[N+d::7])
    c0+=car_car*c_c
    for i in range(n):
      b0[N+d+i*7]+=car_car*b_c[i]
      for j in range(n):
        A0[N+d+i*7,N+d+j*7]+=car_car*A_c[i,j]
  
  # Add in diff constraints for CAR variables which relate adjacent days to each other
  A_c,b_c,c_c=difflogmat(N,order,xx0[N:])
  A0[N:,N:]+=car_car_d*A_c
  b0[N:]+=car_car_d*b_c
  c0+=car_car_d*c_c
  
  # Terms al[i].(I[i]-c[i'].casedata[i'])^2 correspond to CAR error at incidence i (casedata i')
  #       al[i]=( V[ I[i]-c[i'].casedata[i'] ] )^{-1}
  #           ~=(I0[i]^2/(casedata[i'].inc_case))^{-1}
  #             where I0[i] is some approximation to the incidence
  # Also add in normalisation for P(casedata|internal variables xx = incidence, CAR)
  # which means multiplying by prod(xx[N+j]/xx[i])*inc_case^(numcases/2)
  # log(xx[N+j]/xx[i]) = log(xx0[N+j]/xx0[i]) + 2*(xx[N+j]/xx0[N+j]) - (1/2)*(xx[N+j]/xx0[N+j])^2
  #                                           - 2*(xx[i]/xx0[i])     + (1/2)*(xx[i]/xx0[i])^2 + (third order)
  #                   cf -1/2(Ax^2-2bx+c)
  al=np.zeros(N)
  for j in range(len(casedata)):
    if startdate+j not in ignore:
      day=(startdate+j-monday)%7
      i=j-back[day]
      if i>=0:
        al[i]=casedata[j]*inc_case/xx0[i]**2
        A1[i,i]+=al[i]
        A1[i,N+j]-=al[i]*casedata[j]
        A1[N+j,i]-=al[i]*casedata[j]
        A1[N+j,N+j]+=al[i]*casedata[j]**2
        
        A1[N+j,N+j]+=1/xx0[N+j]**2
        A1[i,i]-=1/xx0[i]**2
        b1[N+j]+=2/xx0[N+j]
        b1[i]-=2/xx0[i]
        c1-=2*log(xx0[N+j]/xx0[i])
        numcases+=1
  c1-=numcases*log(inc_case)

  return A0,b0,c0,A1,b1,c1

def getest(enddate=apiday(),prlev=0,eps=1e-3):
  
  N,casedata,onsprev=getextdata(enddate,prlev)
  xx=initialguess(N,onsprev,casedata,back)
  #savevars(N,casedata,back,xx,name="tempinit")
  
  for it in range(20):
    if prlev>=2: print("Iteration",it)
    
    A0,b0,c0,A1,b1,c1=getqform(N,casedata,onsprev,xx)
    
    xx0=xx
    xx=np.maximum(np.linalg.solve(A0+A1,b0+b1),0.01)
    if np.abs(xx/xx0-1).max()<eps: break
  
  #savevars(N,casedata,back,xx,name="England")
  A0,b0,c0,A1,b1,c1=getqform(N,casedata,onsprev,xx)
  A,b,c=A0+A1,b0+b1,c0+c1
  return casedata,xx,A,b,c
  
def getprob(enddate=apiday(),prlev=0,eps=1e-3):
  
  N,casedata,onsprev=getextdata(enddate,prlev)
  xx=initialguess(N,onsprev,casedata,back)
  #savevars(N,casedata,back,xx,name="tempinit")
  
  for it in range(20):
    if prlev>=2: print("Iteration(numerator)",it)
    
    A0,b0,c0,A1,b1,c1=getqform(N,casedata,onsprev,xx)
    
    xx0=xx
    xx=np.maximum(np.linalg.solve(A0+A1,b0+b1),0.01)
    if np.abs(xx/xx0-1).max()<eps: break
  
  num=(-1/2)*np.linalg.slogdet(A0+A1)[1]+(1/2)*(b0+b1)@xx-(1/2)*(c0+c1)
  
  #savevars(N,casedata,back,xx,name="England")

  eta=0.01
  for it in range(20):
    if prlev>=2: print("Iteration(denominator)",it)
    
    A0,b0,c0,A1,b1,c1=getqform(N,casedata,onsprev,xx)
    A0e=A0+eta*np.identity(2*N)
    
    xx0=xx
    xx=np.maximum(np.linalg.solve(A0e,b0),0.01)
    if np.abs(xx/xx0-1).max()<eps: break
  
  denom=(-1/2)*np.linalg.slogdet(A0e)[1]+(1/2)*b0@xx-(1/2)*c0

  return num-denom

# Not done daily CAR-CAR interactions yet
def directeval(xx,casedata,enddate=apiday(),prlev=0):
  enddate=Date(enddate)
  onsprev=getdailyprevalence(maxdate=enddate)
  N=enddate-startdate
  assert len(casedata)<=N

  tx=np.log(xx[:N])
  for o in range(order): tx=tx[:-1]-tx[1:]
  t_II=inc_inc*(tx**2).sum()
  
  # Link incidence to ONSprevalence
  t_IO=0
  for date,targ,sd in onsprev:
    targ/=scale
    sd/=scale
    var=sd**2/inc_ons
    i=date-numk-startdate
    if i>=0 and i+numk<=N:
      t_IO+=((poskern*xx[i:i+numk]).sum()-targ)**2/var
  
  # Add in diff constraints for CAR variables which relate same days of week to each other
  # (d isn't necessarily equal to the day of the week)
  t_CC=0
  for d in range(7):
    n=(N-d+6)//7
    tx=np.log(xx[N+d::7])
    for o in range(order): tx=tx[:-1]-tx[1:]
    t_CC+=car_car*(tx**2).sum()
  
  # Terms al[i].(I[i]-c[i'].casedata[i'])^2 correspond to CAR error at incidence i (casedata i')
  #       al[i]=( V[ I[i]-c[i'].casedata[i'] ] )^{-1}
  #           ~=(I0[i]^2/(casedata[i'].inc_case))^{-1}
  #             where I0[i] is some approximation to the incidence
  t_IC=0
  for j in range(len(casedata)):
    if startdate+j not in ignore:
      day=(startdate+j-monday)%7
      i=j-back[day]
      if i>=0:
        al=casedata[j]*inc_case/xx[i]**2
        t_IC+=al*(xx[N+j]*casedata[j]-xx[i])**2
  
  print("%10.2f %10.2f %10.2f %10.2f       %10.2f"%(t_II,t_IO,t_CC,t_IC,t_II+t_IO+t_CC+t_IC))
  
if 0:
  # Optimising coupling parameters for consistency
  
  def rnd(): return random()*2-1
  seed(42)
  
  while 1:
    inc_ons=exp(-3.2+rnd()*0.5)
    inc_case=exp(0.8+rnd()*1)
    inc_inc=exp(6.2+rnd()*0.8)
    car_car=exp(3.7+rnd()*1)
  
    casedata0,xx0,A0,b0,c0=getest()
    N=xx0.shape[0]//2
  
    now=apiday()
    numcheck=30
    chrange=7
    err=0
    for ch in range(numcheck):
      casedata,xx,A,b,c=getest(now-(ch+1)*chrange,prlev=0)
      i0,i1=N-(ch+2)*chrange,N-(ch+1)*chrange
      #print(xx0[i0:i1])
      #print(xx[i0:i1])
      #print()
      err+=(np.log(xx[i0:i1]/xx0[i0:i1])**2).sum()
      #xx0=xx
    
    print("%12g %12g %12g %12g     %10.6f"%(inc_ons,inc_case,inc_inc,car_car,sqrt(err/(numcheck*chrange))))
    sys.stdout.flush()

if 0:
  # Optimising coupling parameters for consistency - version with natural model
  
  def rnd(): return random()*2-1
  #seed(42)
  
  while 1:
    inc_ons=exp(-3)
    inc_case=exp(3+rnd()*0.2)
    inc_inc=exp(10.9+rnd()*0.2)
    car_car=exp(6+rnd()*0.5)
    car_car_d=exp(2+rnd()*2)
  
    casedata0,xx0,A0,b0,c0=getest(prlev=0)
    N0=xx0.shape[0]//2
    # Interleave to combine incidence and CAR variables so that indices are in date order
    xx0i=np.zeros(2*N0)
    xx0i[0::2]=xx0[:N0]
    xx0i[1::2]=xx0[N0:]
  
    now=apiday()
    numcheck=30
    chrange=7
    dof=LL0=LL1=0
    sys.stdout.flush()
    for ch in range(numcheck):
      casedata,xx,A,b,c=getest(now-(ch+1)*chrange,prlev=0)
      N=N0-(ch+1)*chrange
      assert N==xx.shape[0]//2
      
      # Interleave to combine incidence and CAR variables so that indices are in date order
      xxi=np.zeros(2*N)
      xxi[0::2]=xx[:N]
      xxi[1::2]=xx[N:]
      Ai=np.zeros([2*N,2*N])
      Ai[0::2,0::2]=A[:N,:N]
      Ai[1::2,0::2]=A[N:,:N]
      Ai[0::2,1::2]=A[:N,N:]
      Ai[1::2,1::2]=A[N:,N:]
  
      yyi=xx0i[:2*N]-xxi
      # yy is an observation from the centred normal distribution with precision matrix A
      # but we're only going to look at the bits 2*(N-chrange):2*N, marginalising over the rest.
      # The split is like this:
      # Ai  = [ A_00 A_01 ]
      #       [ A_10 A_11 ]
      # yyi = [ y_0 y_1]
  
      m=2*(N-chrange)
      y_0=yyi[:m]# Not used
      y_1=yyi[m:]
      A_00=Ai[:m,:m]
      A_01=Ai[:m,m:]
      A_10=Ai[m:,:m]
      A_11=Ai[m:,m:]

      # Calculate the marginal precision matrix of the variables indexed from m onwards.
      # Could do C=np.linalg.inv(Ai);B_11=np.linalg.inv(C[m:,m:]), but the following is faster:
      B_11=A_11-A_10@np.linalg.solve(A_00,A_01)
      
      resid=y_1@(B_11@y_1)
      #LL+=(1/2)*np.linalg.slogdet(B_11)[1]-(1/2)*resid
      LL0+=(1/2)*np.linalg.slogdet(B_11)[1]
      LL1+=-(1/2)*resid
      dof+=(1/2)*2*chrange

    #lam=-dof/LL1
    lam=1
    inc_ons*=lam
    inc_case*=lam
    inc_inc*=lam
    car_car*=lam
    LL1*=lam
    LL0+=dof*log(lam)
    print("%12g %12g %12g %12g %12g    %10.6f"%(inc_ons,inc_case,inc_inc,car_car,car_car_d,LL0+LL1))
    sys.stdout.flush()

if 0:
  # Finding overall coupling so as to get a probability distribution over outcomes
  # Deprecated, because now believe we need to fix ONS prevalence coupling externally, not just make self-consistent
  
  casedata0,xx0,A0,b0,c0=getest()
  N0=xx0.shape[0]//2
  # Interleave to combine incidence and CAR variables so that indices are in date order
  xx0i=np.zeros(2*N0)
  xx0i[0::2]=xx0[:N0]
  xx0i[1::2]=xx0[N0:]

  now=apiday()
  numcheck=30
  chrange=7
  tresid=dof=0
  for ch in range(numcheck):
    casedata,xx,A,b,c=getest(now-(ch+1)*chrange,prlev=0)
    N=N0-(ch+1)*chrange
    assert N==xx.shape[0]//2
    # Interleave to combine incidence and CAR variables so that indices are in date order
    xxi=np.zeros(2*N)
    xxi[0::2]=xx[:N]
    xxi[1::2]=xx[N:]
    Ai=np.zeros([2*N,2*N])
    Ai[0::2,0::2]=A[:N,:N]
    Ai[1::2,0::2]=A[N:,:N]
    Ai[0::2,1::2]=A[:N,N:]
    Ai[1::2,1::2]=A[N:,N:]

    yyi=xx0i[:2*N]-xxi
    # yy is an observation from the centred normal distribution with precision matrix A
    # but we're only going to look at the bits 2*(N-chrange):2*N, marginalising over the rest.
    # The split is like this:
    # Ai  = [ A_00 A_01 ]
    #       [ A_10 A_11 ]
    # yyi = [ y_0 y_1]

    m=2*(N-chrange)
    y_0=yyi[:m]# Not used
    y_1=yyi[m:]
    A_00=Ai[:m,:m]
    A_01=Ai[:m,m:]
    A_10=Ai[m:,:m]
    A_11=Ai[m:,m:]

    # Faster version of C=np.linalg.inv(Ai);resid=y_1@(np.linalg.solve(C[m:,m:],y_1)), because don't have to invert the whole matrix:
    resid=y_1@(A_11@y_1)-(y_1@A_10)@np.linalg.solve(A_00,(A_01@y_1))
    print("Resid =",resid/(2*chrange))
    tresid+=resid
    dof+=2*chrange
  lam=tresid/dof
  print("Overall residual factor =",lam)
    
if 0:
  seed(42)
  np.random.seed(42)
  casedata,xx0,A,b,c=getest(prlev=2)
  N=xx0.shape[0]//2
  # getcaseoutliers(casedata,N)
  conf=0.95
  C=np.linalg.inv(A)
  nsamp=10000
  l=mvn.rvs(mean=xx0,cov=C,size=nsamp)
  l[:,:N]=np.maximum(l[:,:N],0.01)
  l[:,N:]=np.maximum(l[:,N:],1)
  np.ndarray.sort(l,axis=0)
  i0=int(nsamp*(1-conf)/2)
  i1=int(nsamp*(1+conf)/2)
  low=l[i0,:]
  high=l[i1,:]
  savevars(N,casedata,back,xx0,low=low,high=high,name="England")
  poi
  
if 0:
  inc_ons=exp(5.0)
  inc_case=exp(20)
  inc_inc=exp(8.9)
  car_car=exp(-0.5)
  car_car_d=exp(3)
  
  inc_ons=exp(-3+rnd()*0)
  inc_case=exp(9+rnd()*1.5)
  inc_inc=exp(6.5+rnd()*1)
  car_car=exp(8+rnd()*2)
  car_car_d=exp(0.+rnd()*2)
  
  inc_ons=exp(-3+rnd()*0)
  inc_case=exp(0+rnd()*0)
  inc_inc=exp(8.5)
  car_car=exp(-0.5)
  car_car_d=exp(-0.5)
    
  casedata,xx0,A,b,c=getest(prlev=2)
  N=xx0.shape[0]//2
  # getcaseoutliers(casedata,N)
  C=np.linalg.inv(A)
  nsamp=20
  l=mvn.rvs(mean=xx0,cov=C,size=nsamp)
  l[:,:N]=np.maximum(l[:,:N],0.01)
  l[:,N:]=np.maximum(l[:,N:],1)
  print(xx0[N-14:N],xx0[-14:],end="")
  directeval(xx0,casedata)
  print()
  for i in range(nsamp):
    print(l[i][N-14:N],l[i][-14:],end="")
    directeval(l[i],casedata)
  print()
  print(xx0[:N-5]/xx0[N+5:]/casedata[5:])
  poi

if 1:
  while 1:
    inc_ons=exp(-3+rnd()*0)
    inc_case=exp(0+rnd()*0)
    inc_inc=exp(8.5+rnd()*1)
    car_car=exp(-0.5+rnd()*2)
    car_car_d=exp(-0.5+rnd()*2)
    
    inc_ons=exp(-3+rnd()*0)
    inc_case=exp(0+rnd()*0)
    inc_inc=exp(8.5+rnd()*0)
    car_car=exp(-0.5+rnd()*0)
    car_car_d=exp(-0.5+rnd()*0)
    
    inc_ons=exp(5.0+rnd()*0.)
    inc_case=exp(10+rnd()*0)
    inc_inc=exp(8.9+rnd()*0.0)
    car_car=exp(-0.5+rnd()*0.0)
    car_car_d=exp(3+rnd()*0.0)
    
    inc_ons=exp(3.5+rnd()*1.5)# 13.5
    inc_case=exp(4+rnd()*2)
    inc_inc=exp(7.5+rnd()*2)
    car_car=exp(5.5+rnd()*2)
    car_car_d=exp(2+rnd()*2)
    #2.54905e+06      38.2607      3753.01      1164.38      3.61769    11721398.574450

    inc_ons=exp(13.5+rnd())# 13.5
    inc_case=exp(4+rnd()*0)
    inc_inc=exp(7.5+rnd()*0)
    car_car=exp(5.5+rnd()*0)
    car_car_d=exp(2+rnd()*0)
    
    LL=getprob(prlev=0)
    print("%12g %12g %12g %12g %12g    %10.6f"%(inc_ons,inc_case,inc_inc,car_car,car_car_d,LL))
    sys.stdout.flush()
    
