# Try out simulation of areas where a significant proportion of the population caught it,
# and there has been an antibody survey that gives some kind of decent prevalence
# estimate.
# This way, should be able to get (a weak) lower bound on disease-induced HIT under Gomes modification.

import csv,sys,getdata,requests
if sys.version_info[0]<3: print("Requires Python 3");sys.exit(1)
import numpy as np
from scipy.special import gammainc
from scipy.stats import gamma as gammadist

import scipy
import scipy.special
def fact(x): return scipy.special.gamma(x+1)
# If n is negative then r must be a non-negative integer in this version (which is true here)
def bin(n,r):
  if n>=-0.6:
    return fact(n)/(fact(r)*fact(n-r))
  else:
    return (-1)**r*bin(r-n-1,r)

from math import log,exp

location="London"

# subdivisions of a day (not used yet)
# subd=10

# Return a negative binomial distribution of mean mu, variance var, clipped to [mi,mx]
def getNB(mu,var,mi,mx):
  p=1-mu/var
  assert p>0 and p<1
  r=mu*(1-p)/p
  dist=np.zeros(mx+1)
  for k in range(mi,mx+1): dist[k]=bin(k+r-1,k)*p**k
  s=dist.sum();dist/=s
  return dist

def getGamma(k,maxsbins):
  # Space out bins according to shape k+1 so that each bin represents an equal
  # (incomplete) expectation of X in the shape k distribution. Other bin choices are
  # possible, but this is a good one from the point of view of achieving more accuracy
  # using a smaller number of bins.
  l=gammadist.ppf([i/maxsbins for i in range(maxsbins)],k+1)

  # Calculate:
  #   m0[i]   = P[X < l_i]
  #   m1[i]   = E[X; X < l_i]
  #   susc[i] = E[X | l_i <= X < l_{i+1}], the representative susceptibility for bin i
  #   q[i]    = P(l_i <= X < l_{i+1})
  m0=np.append(gammainc(k,l),1)
  m1=np.append(gammainc(k+1,l),1)
  susc=np.array([(m1[i+1]-m1[i])/(m0[i+1]-m0[i]) for i in range(maxsbins)])
  q=np.array([m0[i+1]-m0[i] for i in range(maxsbins)])
  
  return susc,q

def getConst():
  return np.array([1.]),np.array([1.])

def getqConst():
  return np.array([0.99999,1.00001]),np.array([0.5,0.5])

# YYYY-MM-DD -> day number
def datetoday(s):
  mm=[0,31,59,90,120,151,181,212,243,273,304,334]
  y=int(s[:4])
  m=int(s[5:7])
  d=int(s[8:10])
  return (y-1970)*365+(y-1969)//4+mm[m-1]+(m>=3 and (y&3)==0)+d-1

def daytodate(n):
  mm=[31,28,31,30,31,30,31,31,30,31,30,31]
  y=1970+(n//1461)*4;n%=1461
  if n>=365: n-=365;y+=1
  if n>=365: n-=365;y+=1
  if n>=366: n-=366;y+=1
  m=0
  while 1:
    o=mm[m]+(m==1 and (y&3)==0)
    if n<o: break
    m+=1;n-=o
  return "%4d-%02d-%02d"%(y,m+1,n+1)

# Minimum, maximum time from being infected to being infectious; ditto reporting cases
mini=2;maxi=14
minr=0;maxr=100

# Estimated number of days from infection for antibodies to show up, plus time to report this in a survey.
# This actually makes very little difference in the cases of NYC and London, because the antibody surveys were very
# late relative to the peak - i.e., the rate of infection was almost negligible at the time of survey.
ABtime=21

if location=="NYC":
  # NYC case count from https://www.worldometers.info/coronavirus/usa/new-york/
  repcasestart=datetoday("2020-03-13")
  repcases=np.array([93, 107, 212, 235, 742, 1342, 2341, 3052, 1993, 5440, 5123, 5516, 6674, 6097, 7380, 7250, 7413, 6785, 8823, 8104, 9353, 10628, 11506, 8477, 9135, 10714, 9000, 10533, 11045, 8969, 8458, 6419, 7614, 11661, 7636, 7753, 7090, 6174, 4879, 4461, 5713, 6313, 8864, 10868, 5678, 4013, 3446, 4708, 4681, 4383, 3991, 4670, 3491, 2765, 3352, 3930, 3284, 2704, 1997, 1745, 1504, 2193, 2248, 2920, 2083, 1748, 1419, 1364, 1619, 2108, 1579, 1720, 1537, 1301],dtype=float)
  pop=8.4e6
  start  = datetoday("2020-03-01")
  change = datetoday("2020-03-19")
  end    = datetoday("2020-06-15")
  # https://www.nbcnewyork.com/news/local/cuomo-outlines-reopening-roadmap-for-new-york-as-daily-deaths-hit-lowest-level-in-weeks/2390949/
  survey=0.247
  surveyday=max(datetoday("2020-04-27")-ABtime,start)
  
elif location=="London":
  if 0:
    # From https://en.wikipedia.org/wiki/COVID-19_pandemic_in_London#Data
    repcasestart=datetoday("2020-03-11")
    repcases=np.array([13, 32, 31, 146, 94, 73, 141, 332, 268, 367, 377, 224, 244, 439, 375, 672, 718, 662, 658, 564, 600, 1220, 950, 956, 517, 1214, 658, 742, 977, 916, 740, 710, 758, 521, 472, 479, 560, 704, 689, 453, 297, 418, 280, 415, 296, 278, 267, 225, 146, 111, 207, 180, 223, 128, 160, 252, 117, 142, 124, 111, 95, 62, 89, 44, 89, 80, 166, 81, 47, 42, 50, 49, 55, 59, 38, 38, 22, 27, 16, 30, 25, 30, 21, 32, 21, 21, 36],dtype=float)
  else:
    # From https://coronavirus.data.gov.uk/
    url="https://coronavirus.data.gov.uk/downloads/csv/coronavirus-cases_latest.csv"
    with requests.get(url) as resp:
      if resp.status_code!=200: raise ConnectionError("Couldn't load "+url)
      r=csv.reader(resp.content.decode('utf-8').splitlines())
      # Fields = ['Area name', 'Area code', 'Area type', 'Specimen date', 'Daily lab-confirmed cases', 'Previously reported daily cases', 'Change in daily cases', 'Cumulative lab-confirmed cases', 'Previously reported cumulative cases', 'Change in cumulative cases', 'Cumulative lab-confirmed cases rate']
    d={};repcasestart=1e30;mx=-1e30
    for x in r:
      if x[0]=='London':
        i=datetoday(x[3])
        d[i]=int(float(x[4]))
        if i<repcasestart: repcasestart=i
        if i>mx: mx=i
    l=[]
    for j in range(i,mx+1): l.append(d.get(j,0))
    repcases=np.array(l,dtype=float)
  pop=8.982e6
  start  = datetoday("2020-03-01")
  change0 = datetoday("2020-03-23")
  end    = datetoday("2020-06-15")
  # https://www.itv.com/news/2020-05-21/health-secretary-matt-hancock-government-daily-coronavirus-press-conference/
  survey=0.17
  surveyday=max(datetoday("2020-05-21")-ABtime,start)
else:
  raise LookupError("Unrecognised location: %s"%location)

#o=datetoday(repcasestart)
#for (i,c) in enumerate(repcases):
#  print(daytodate(o+i),"%6d"%c)

# Parameters to be stochastically varied in later version
# Using guesstimate numbers to play with to start off with.
change=change0-2
k=2 # connectivity dispersion parameter; lower = more dispersed
R0a=3 # pre-lockdown
R0b=0.65 # post-lockdown
infectmean=4
infectvar=7
if location=="London": reportmean=8;reportvar=80
else: reportmean=25;reportvar=250

nbins=50
infectdist=getNB(infectmean,infectvar,mini,maxi)
reportdist=getNB(reportmean,reportvar,minr,maxr)
#conn,pconn=getConst()
#conn,pconn=getqConst()
conn,pconn=getGamma(k,nbins)

def getnewcases(conn,pconn,start,change,end,initial,pop,R0a,R0b,infectdist):
  # new[len(new)-1-i] = simulated number of people who became (newly) infected exactly i days ago (vectorised over suscep)
  if initial>pop: initial=pop
  pre=5;new=[initial/pre*pconn]*pre# simplification pro tem for the pre-historical segment
  suscep=pop*pconn-sum(new)
  for date in range(start,end+1):
    t=0
    for j in range(mini,maxi+1):
      if j>len(new): break
      t+=(infectdist[j]*new[len(new)-j]*conn).sum()# This factor of conn reflects different infectivities
    R=(R0a if date<change else R0b)
    new.append(np.minimum(R*t*conn*suscep/pop,suscep))# This factor of conn reflects different susceptibilities
    suscep-=new[-1]
  #print(initial,(np.array([x.sum() for x in new[pre:]])[:surveyday-start]).sum())
  return np.array([x.sum() for x in new[pre:]])

# Find the initial count that gives rise to the same total prevalence as in the survey
target=survey*pop
initial0=10000
n0=getnewcases(conn,pconn,start,change,end,initial0,pop,R0a,R0b,infectdist)[:surveyday-start].sum()# alter: add initial
if n0<target:
  while 1:
    initial1=initial0*2
    n1=getnewcases(conn,pconn,start,change,end,initial1,pop,R0a,R0b,infectdist)[:surveyday-start].sum()
    if n1>=target: break
    initial0=initial1;n0=n1
else:
  while 1:
    initial1=initial0;n1=n0
    initial0=initial1/2
    n0=getnewcases(conn,pconn,start,change,end,initial0,pop,R0a,R0b,infectdist)[:surveyday-start].sum()
    if n0<target: break
ff=0.1
while abs(n1/n0-1)>1e-6:
  initial=(target-n0)/(n1-n0)*(initial1-initial0)+initial0
  initial0=initial-ff*(initial-initial0)
  initial1=initial+ff*(initial1-initial)
  n0=getnewcases(conn,pconn,start,change,end,initial0,pop,R0a,R0b,infectdist)[:surveyday-start].sum()
  n1=getnewcases(conn,pconn,start,change,end,initial1,pop,R0a,R0b,infectdist)[:surveyday-start].sum()
initial=(initial0+initial1)/2
print("Initial infections on",daytodate(start),"=",initial,file=sys.stderr)

# new[i]     = predicted number of infections on day start+i
# predrep[i] = predicted number of reported cases on day start+i

new=getnewcases(conn,pconn,start,change,end,initial,pop,R0a,R0b,infectdist)

predrep=[]
for i in range(len(new)):
  t=0
  for j in range(minr,maxr+1):
    if j>i: break
    t+=reportdist[j]*new[i-j]
  predrep.append(t)
predrep=np.array(predrep)

repcaseend=repcasestart+len(repcases)-1
# Intersect [start,end] (represented in predrep[])
# and       [repcasestart,repcaseend] (represented in repcases)
d0=max(start,repcasestart)
d1=min(end,repcaseend)+1
psum=sum(predrep[d0-start:d1-start])
rsum=sum(repcases[d0-repcasestart:d1-repcasestart])

# Scale reported cases to match estimated number of actual cases
repratio=rsum/psum
repcases/=repratio

for day in range(start,end+1):
  print("%4d  %s  %9.5f  %9.5f"%(day-start,daytodate(day),new[day-start]/pop*100,predrep[day-start]/pop*100),end="")
  if day>=repcasestart and day<=repcaseend: print("  %9.5f"%(repcases[day-repcasestart]/pop*100))
  else: print("          -")

# Can't deduce HIT in general because post-lockdown, (high R0, low HIT) looks like (low R0, high HIT).
# But we _can_ probably eliminate some of the madly low HITs because they would require an extreme variation in effective R

