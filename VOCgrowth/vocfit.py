from stuff import *
import sys,re,argparse,pickle
from scipy.optimize import minimize
from math import log,exp,sqrt,sin,pi
import numpy as np
from subprocess import Popen,PIPE

# (Make it auto download files?)
# Get ltla.csv from https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv
# Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv
# COG-UK data from https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv
# SGTF   data from Fig.16 Tech Briefing 12: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/988608/Variants_of_Concern_Technical_Briefing_12_Data_England.xlsx

apicases=loadcsv("ltla.csv")
ltlaengdata=loadcsv("Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv")
ltlaukdata=loadcsv("Local_Authority_District_to_Country_(April_2019)_Lookup_in_the_United_Kingdom.csv")
ltla2ltla=dict(zip(ltlaukdata['LAD19CD'],ltlaukdata['LAD19CD']))
ltla2uk=dict((ltla,"UK") for ltla in ltlaukdata['LAD19CD'])
ltla2country=dict(zip(ltlaukdata['LAD19CD'],ltlaukdata['CTRY19NM']))
ltla2region=dict(ltla2country,**dict(zip(ltlaengdata['LAD19CD'],ltlaengdata['RGN19NM'])))
def coglab2uk(x): return "UK"
def coglab2country(x): return x.split('/')[0].replace('_',' ')
def coglab2coglab(x): return x
def sgtf2region(x):
  if x=='Yorkshire and Humber': return 'Yorkshire and The Humber'
  return x
def sgtf2country(x): return 'England'
def includeltla(ltla,ltlaset):
  if ltlaset=="London":
    return ltla2region[ltla]=='London'
  elif ltlaset=="Bolton":
    return ltla=='E08000001'
  elif ltlaset=="All":
    return True
  else:
    raise RuntimeError("Unrecognised ltla set "+ltlaset)

parser=argparse.ArgumentParser()
parser.add_argument('-l', '--load-options', help='Load options from a file')
parser.add_argument('-s', '--save-options', help='Save options to a file')
parser.add_argument('-g', '--graph-filename', help='Stem of graph filenames')
args=parser.parse_args()

### Model ###
#
# Known:
# n_i      = number of confirmed cases on day i by specimen date (slightly adjusted for weekday)
# p        = case ascertainment rate (chance of seeing a case)
# g_{-1}(r),v_{-1}(r) = emperical growth rate (and its variance) in the two weeks up to 2021-04-10 as a function of region, r
# r_j, s_j = Variant counts of non-B.1.617.2, B.1.617.2 in j^th week
# I_j      = set of days (week) corresponding to VOC counts r_j, s_j
# Assume chance of sequencing a case is a totally free parameter, and optimise over it
#
# Unknown:
# h   = daily growth advantage of B.1.617.2 over other variants
# X_n = Fourier coefficients controlling the growth of non-B.1.617.2
# A_0 = initial count of non-B.1.617.2
# B_0 = initial count of B.1.617.2
#
# Likelihood:
# A_{i+1}=e^{g_i}A_i
# B_{i+1}=e^{g_i+h}B_i
# n_i ~ Po(p(A_i+B_i))
# r_j ~ Po(q_j.A_{I_j})  (A_{I_j} means sum_{i in A_j}A_i)
# s_j ~ Po(q_j.B_{I_j})
# [q_j ~ eps.e^{-eps.q}dq (a wide exponential distribution) - not currently implemented]
# q_j is optimised out, so q_j=(r_j+s_j)/(A_{I_j}+B_{I_j}), or effectively A_{I_j}/(A_{I_j}+B_{I_j}) ~ Beta(r_j+1,s_j+1)
# g_0 ~ N(g_{-1},v_{-1})
# X_n ~ N(0,1)
# N=ndays-1
# g_i = g_0 + bmscale*sqrt(N)*(i/N*X_0 + sqrt(2)/pi*sum_n e^{-(n*bmsig/N)^2/2}sin(n*pi*i/N)*X_n/n)
# h ~ N(0,tau^2)
#
### End Model ###


### Options ###

#source="Sanger"
source="COG-UK"
#source="SGTF"

# Can choose location size from "LTLA", "region", "country", "UK"
# Sanger works with LTLA, region, country
# COG-UK works with country, UK
# SGTF works with region, country
locationsize="UK"

ltlaexclude=set()
# ltlaexclude=set(['E08000001'])# This would exclude Bolton
ltlaset="All"
#ltlaset="London"
#ltlaset="Bolton"

mgt=5# Mean generation time in days

# Earliest day to use case data
minday=datetoday('2021-04-01')# Inclusive

# Earliest day to use VOC count data, given as end-of-week. Will be rounded up to match same day of week as lastweek.
#firstweek=minday+6
firstweek=datetoday('2021-04-17')

nif1=1  # Non-independence factor for cases (less than 1 means downweight this information)
nif2=1  # Non-independence factor for VOC counts (ditto)
isd=1   # Inverse sd for prior on transmission advantage (as growth rate per day). 0 means uniform prior. 1 is very weak.

# Prior linking initial daily growth rate to estimate from pre-B.1.617.2 era
sig0=0.004

# Timescale in days over which growth rate can change significantly
# (lower = more wiggles)
bmsig=30

# Lengthscale for filtered Brownian motion
# (higher = greater amplitude for the wiggles)
bmscale=0.02

# Case ascertainment rate
asc=0.4

# Discard this many cases at the end of the list of cases by specimen day
discardcasedays=2

# Discard this many days of the latest COG data
discardcogdays=2

# Collect together all locations without positive entries into one combined "Other" location
# (Makes little difference in practice)
bundleremainder=True

minopts={"maxiter":10000,"eps":1e-5}

mode="local growth rates"
#mode="global growth rate"
#mode="fixed growth rate",0.1

voclen=(1 if source=="COG-UK" else 7)

### End options ###

opts={
  "Source": source,
  "Location size": locationsize,
  "LTLA set": ltlaset,
  "LTLA exclude": list(ltlaexclude),
  "Generation time (days)": mgt,
  "Earliest day for case data": daytodate(minday),
  "Earliest week (using end of week date) to use VOC count data": daytodate(firstweek),
  "Optimisation mode": mode,
  "nif1": nif1,
  "nif2": nif2,
  "Inverse sd for prior on growth": isd,
  "Sigma (prior on daily growth rate change)": sig0,
  "Timescale of growth rate change (days)": bmsig,
  "Lengthscale for filtered Brownian motion": bmscale,
  "Case ascertainment rate": asc,
  "Number of days of case data to discard": discardcasedays,
  "Number of days of COG-UK data to discard": discardcogdays,
  "Bundle remainder": bundleremainder,
  "Minimiser options": minopts,
  "Length of time period over which VOC counts are given (days)": voclen
}

if args.save_options!=None:
  with open(args.save_options,'w') as fp: json.dump(opts,fp,indent=2)

if args.load_options!=None:
  with open(args.load_options,'r') as fp: lopts=json.load(fp)
  for x in lopts: opts[x]=lopts[x]
  
source=opts["Source"]
locationsize=opts["Location size"]
ltlaset=opts["LTLA set"]
ltlaexclude=set(opts["LTLA exclude"])
mgt=opts["Generation time (days)"]
minday=datetoday(opts["Earliest day for case data"])
firstweek=datetoday(opts["Earliest week (using end of week date) to use VOC count data"])
mode=opts["Optimisation mode"]
nif1=opts["nif1"]
nif2=opts["nif2"]
isd=opts["Inverse sd for prior on growth"]
sig0=opts["Sigma (prior on daily growth rate change)"]
bmsig=opts["Timescale of growth rate change (days)"]
bmscale=opts["Lengthscale for filtered Brownian motion"]
asc=opts["Case ascertainment rate"]
discardcasedays=opts["Number of days of case data to discard"]
discardcogdays=opts["Number of days of COG-UK data to discard"]
bundleremainder=opts["Bundle remainder"]
minopts=opts["Minimiser options"]
voclen=opts["Length of time period over which VOC counts are given (days)"]

print("Options:")
print()
for x in sorted(list(opts)): print("%s:"%x,opts[x])
print()
sys.stdout.flush()

np.set_printoptions(precision=3,linewidth=150)

if ltlaset=="All":
  if source=="COG-UK": areacovered="the UK"
  else: areacovered="England"
else:
  areacovered=ltlaset

maxday=datetoday(max(apicases['date']))-discardcasedays# Inclusive
ndays=maxday-minday+1
firstweek=max(firstweek,minday+voclen-1)

if source=="Sanger":
  fullsource="Wellcome Sanger Institute"
  assert voclen==7
  sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
  
  lastweek=datetoday(max(sanger['WeekEndDate']));assert maxday>=lastweek
  nweeks=(lastweek-firstweek)//voclen+1
  # Sanger week number is nweeks-1-(lastweek-day)//voclen

  if locationsize=="LTLA":
    reduceltla=ltla2ltla
  elif locationsize=="region":
    reduceltla=ltla2region
  elif locationsize=="country":
    reduceltla=ltla2country
  else:
    raise RuntimeError("Incompatible source, locationsize combination: "+source+", "+locationsize)
  
  # Get Sanger (variant) data into a suitable form
  vocnum={}
  for (date,ltla,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
    if ltla in ltlaexclude or not includeltla(ltla,ltlaset): continue
    day=datetoday(date)
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks:
      place=reduceltla[ltla]
      if place not in vocnum: vocnum[place]=np.zeros([nweeks,2],dtype=int)
      if var=="B.1.617.2": vocnum[place][week][1]+=n
      else: vocnum[place][week][0]+=n
elif source=="COG-UK":
  fullsource="COG-UK"
  cog=loadcsv("cog_metadata.csv")
  lastweek=datetoday(max(cog['sample_date']))-discardcogdays;assert maxday>=lastweek
  nweeks=(lastweek-firstweek)//voclen+1
  # Week number is nweeks-1-(lastweek-day)//voclen
  
  if  locationsize=="country":
    reduceltla=ltla2country
    reducecog=coglab2country
  elif locationsize=="UK":
    reduceltla=ltla2uk
    reducecog=coglab2uk
  else:
    raise RuntimeError("Incompatible source, locationsize combination: "+source+", "+locationsize)
  
  # Get COG-UK (variant) data into a suitable form
  vocnum={}
  for (date,seqname,var) in zip(cog['sample_date'],cog['sequence_name'],cog['lineage']):
    day=datetoday(date)
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks:
      r=re.match("[^0-9-]*[0-9-]",seqname)
      coglab=seqname[:r.end()-1]
      place=reducecog(coglab)
      if place not in vocnum: vocnum[place]=np.zeros([nweeks,2],dtype=int)
      if var=="B.1.617.2": vocnum[place][week][1]+=1
      else: vocnum[place][week][0]+=1
elif source=="SGTF":
  fullsource="SGTF data from PHE Technical briefing 13"
  assert voclen==7
  sgtf=loadcsv("TechBriefing13Fig19.csv")
  lastweek=max(datetoday(x) for x in sgtf['week'])+6# Convert w/c to w/e convention
  assert maxday>=lastweek
  nweeks=(lastweek-firstweek)//voclen+1
  # Week number is nweeks-1-(lastweek-day)//voclen

  if locationsize=="region":
    reduceltla=ltla2region
    reducesgtf=sgtf2region
  elif  locationsize=="country":
    reduceltla=ltla2country
    reducesgtf=sgtf2country
  else:
    raise RuntimeError("Incompatible source, locationsize combination: "+source+", "+locationsize)
  
  # Get SGTF data into a suitable form
  vocnum={}
  for (date,region,var,n) in zip(sgtf['week'],sgtf['Region'],sgtf['sgtf'],sgtf['n']):
    day=datetoday(date)+6# Convert from w/c to w/e convention
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks:
      place=reducesgtf(region)
      if place not in vocnum: vocnum[place]=np.zeros([nweeks,2],dtype=int)
      vocnum[place][week][int("SGTF" not in var)]+=n
  # Adjust for non-B.1.617.2 S gene positives, based on the assumption that these are in a non-location-dependent proportion to the number of B.1.1.7
  # This is likely to be dodgy
  # From COG-UK: B117  Others (not B.1.617.2)
  # 2021-03-11  18295     341   1.86%
  # 2021-03-18  16900     308   1.82%
  # 2021-03-25  14920     231   1.55%
  # 2021-04-01  10041     243   2.42%
  # 2021-04-08   7514     240   3.19%
  # 2021-04-15   7150     348   4.87%
  # 2021-04-22   6623     367   5.54%
  # 2021-04-29   5029     205   4.08%
  # 2021-05-06   4383     217   4.95%
  date0,date1=datetoday('2021-03-11'),datetoday('2021-04-15')
  for place in vocnum:
    for week in range(nweeks):
      day=lastweek-voclen*(nweeks-1-week)
      assert day>=date0
      if day<date1: f=(day-date0)/(date1-date0)*0.03+0.02
      else: f=0.05
      vocnum[place][week][1]=max(vocnum[place][week][1]-int(f*vocnum[place][week][0]),0)
else:
  raise RuntimeError("Unrecognised source: "+source)

# Simple weekday adjustment by dividing by the average count for that day of the week.
# Use a relatively stable period (inclusive) over which to take the weekday averages.
weekadjdates=[datetoday('2021-04-03'),datetoday('2021-05-14')]
weekadj=np.zeros(7)
for (date,n) in zip(apicases['date'],apicases['newCasesBySpecimenDate']):
  day=datetoday(date)
  if day>=weekadjdates[0] and day<=weekadjdates[1]: weekadj[day%7]+=n
weekadjp=weekadj*7/sum(weekadj)

# Get case data into a suitable form
preweek=minday+9# Gather pre-variant counts in two one-week periods up to this date
precases0={}
cases={}
for (ltla,date,n) in zip(apicases['areaCode'],apicases['date'],apicases['newCasesBySpecimenDate']):
  if ltla not in reduceltla or ltla in ltlaexclude or not includeltla(ltla,ltlaset): continue
  day=datetoday(date)
  d=day-minday
  place=reduceltla[ltla]
  if place not in vocnum: continue
  if place not in precases0: precases0[place]=np.zeros(2,dtype=int)
  if place not in cases: cases[place]=np.zeros(ndays)
  for i in range(2):
    if day>preweek-7*(2-i) and day<=preweek-7*(1-i):
      precases0[place][i]+=n
  if d>=0 and d<ndays:
    cases[place][d]+=n/weekadjp[day%7]
places=sorted(list(cases))

# Restrict to places for which there is at least some of each variant, and bundle the remaining locations together as "Other"
okplaces=set([place for place in places if vocnum[place][:,0].sum()>0 and vocnum[place][:,1].sum()>0])
if bundleremainder:
  otherplaces=set(places).difference(okplaces)
  othervocnum=sum((vocnum[place] for place in otherplaces),np.zeros([nweeks,2],dtype=int))
  othercases=sum((cases[place] for place in otherplaces),np.zeros(ndays,dtype=int))
  if othervocnum[:,0].sum()>0 and othervocnum[:,1].sum()>0:
    okplaces.add("Other")
    vocnum["Other"]=othervocnum
    cases["Other"]=othercases
places=list(okplaces)
places.sort()# Alphabetical order

# Work out pre-B.1.617.2 case counts, amalgamated to at least region level
precases={}
for place in precases0:
  if bundleremainder and place in otherplaces:
    dest="Other"
  elif locationsize=="LTLA": dest=ltla2region[place]
  else: dest=place
  if dest not in precases: precases[dest]=np.zeros(2,dtype=int)
  precases[dest]+=precases0[place]
def prereduce(place):
  if place!="Other" and locationsize=="LTLA": return ltla2region[place]
  else: return place

# Convert daily growth rate & uncertainty into R-number-based description
def Rdesc(h0,dh):
  z=1.96
  (Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-z*dh,h0,h0+z*dh]]
  return "%.0f%% (%.0f%% - %.0f%%)"%(T,Tmin,Tmax)

print("Estimating transmission advantage using variant counts only (not case counts)")
print("=============================================================================")
print()

from scipy.special import betaln
from scipy.integrate import quad
from scipy import inf
def crossratiosubdivide(matgen):
  tot=np.zeros([2,2],dtype=int)
  ndiv=20
  hmin=0;hmax=0.3# Range of possible daily growth advantages
  logp=np.zeros(ndiv)
  L0=L1=0
  for M in matgen:
    tot+=M
    if (M>0).all():
      c=1/((1/M.flatten()).sum())
      T=M[0,0]*M[1,1]/(M[0,1]*M[1,0])
      L0+=c*log(T);L1+=c
      for i in range(ndiv):
        x=(hmin+(i+.5)/ndiv*(hmax-hmin))*voclen# Convert to weekly growth rate
        a,b,c,d=M[0,0],M[0,1],M[1,0],M[1,1]
        l0=d*x-(betaln(a,b)+betaln(c,d))
        # Faff around finding maximum to avoid underflow in integral
        e=exp(x)
        X=b+d-1;Y=a+b;Z=c+d
        A=Y+Z-X
        B=Y/e+Z-X*(1+1/e)
        C=-X/e
        z=(-B+sqrt(B**2-4*A*C))/(2*A)
        l1=(b+d-1)*log(z) - (a+b)*log(1+z) - (c+d)*log(1+e*z)
        res=quad(lambda z: exp( (b+d-1)*log(z) - (a+b)*log(1+z) - (c+d)*log(1+e*z) - l1 ), 0, inf)
        logp[i]+=log(res[0])+l0+l1
  if (tot==0).any():
    print("Can't estimate best transmission factor because VOC count matrix has 1 or more zero entries");return
  g=log(tot[0,0]*tot[1,1]/(tot[0,1]*tot[1,0]))/voclen
  dg=sqrt((1/tot.flatten()).sum())/voclen
  print("Overall cross ratio:",Rdesc(g,dg),tot.flatten())
  print("Inverse variance weighting method using log(CR):",Rdesc(L0/L1/voclen,1/sqrt(L1)/voclen))
  i=np.argmax(logp)
  if i==0 or i==ndiv-1:
    print("Can't properly estimate best transmission factor or confidence interval because the maximum is at the end")
    imax=i
    c=0.1
  else:
    b=(logp[i+1]-logp[i-1])/2
    c=2*logp[i]-(logp[i+1]+logp[i-1])
    imax=i+b/c
  irange=1/sqrt(c)
  g0=(hmin+(hmax-hmin)*(imax+.5)/ndiv)
  dg=(hmax-hmin)*irange/ndiv
  print("Likelihood method using log(CR):",Rdesc(g0,dg))
  print()

if voclen>=7:
  for w in range(nweeks-1):
    day0=lastweek-(nweeks-w)*voclen+1
    print(daytodate(day0),"-",daytodate(day0+2*voclen-1))
    crossratiosubdivide(vocnum[place][w:w+2] for place in places)
print("All week pairs:")
crossratiosubdivide(vocnum[place][w:w+2] for place in places for w in range(nweeks-1))

print("Estimating transmission advantage using variant counts together with case counts")
print("================================================================================")
print()

# L=ndays-1
# i is time from start, in days
# t=i/L
# growth[i] = bmscale*sqrt(L)*(t*X_0 + sum_{n=1}^{N-1} sqrt(2)/pi*exp(-(n*bmsig/L)^2/2)*sin(n*pi*t)/n*X_n)
# where X_n ~ N(0,1),  n=0,...,N; N=ceil(4*L/bmsig), say

bmL=ndays-1
bmN=int(3*bmL/bmsig+1)
bmsin=[sin(r*pi/bmL) for r in range(2*bmL)]
bmweight=[0]+[sqrt(2)/pi*exp(-(n*bmsig/bmL)**2/2)/n for n in range(1,bmN)]

# Need to scale the variables being optimised over to keep SLSQP happy
condition=np.array([50,50,1000,1000]+[1.]*bmN)

# bmN+4 parameters to be optimised:
# 0: a0
# 1: b0
# 2: h
# 3: g0
# 4 ... 4+bmN-1 : X_0, ..., X_{bmN-1}

def expand(xx):
  (a0,b0)=xx[:2]
  AA=[exp(a0)];BB=[exp(b0)]
  h=xx[2]
  g0=xx[3]
  a=a0;b=b0
  w=bmscale*sqrt(bmL)
  GG=[]
  for i in range(ndays-1):
    t=i/bmL
    gu=t*xx[4]
    for n in range(1,bmN):
      gu+=bmweight[n]*bmsin[(i*n)%(2*bmL)]*xx[4+n]
    g=g0+w*gu
    a+=g;b+=g+h;GG.append(g)
    AA.append(exp(a))
    BB.append(exp(b))
  return AA,BB,GG

# Return negative log likelihood
def NLL(xx_conditioned,lcases,lvocnum,sig0,p,lprecases):
  xx=xx_conditioned/condition
  a,b=lprecases[0]+.5,lprecases[1]+.5
  g0=log(b/a)/7
  v0=(1/a+1/b)/49+sig0**2
  tot=-(xx[3]-g0)**2/(2*v0)
  AA,BB,GG=expand(xx)
  # Component of likelihood due to number of confirmed cases seen
  for i in range(ndays):
    lam=p*(AA[i]+BB[i])
    # max with -10000 because the expression is unbounded below which can cause a problem for SLSQP
    tot+=max((lcases[i]-lam+lcases[i]*log(lam))*nif1,-10000)
  # Term to regulate change in growth rate
  for i in range(bmN):
    tot+=-xx[4+i]**2/2
  # Term to align the variant numbers with VOC count data
  for w in range(nweeks):
    endweek=lastweek-(nweeks-1-w)*voclen-minday
    A=sum(AA[endweek-(voclen-1):endweek+1])
    B=sum(BB[endweek-(voclen-1):endweek+1])
    tot+=(lvocnum[w][0]*log(A/(A+B))+lvocnum[w][1]*log(B/(A+B)))*nif2
  # Prior on h
  tot+=-(xx[2]*isd)**2/2
  return -tot

def optimiseplace(place,hint=np.zeros(bmN+4),fixedh=None,statphase=False):
  xx=np.copy(hint)
  # bounds[2][0]=0 prejudges B.1.617.2 as being at least as transmissible as B.1.1.7. This helps SLSQP not get stuck in some cases
  # though would need to relax this constraint if dealing with other variants where it might not be true.
  bounds=[(-10,20),(-10,20),(0,1),(-1,1)]+[(-10,10)]*bmN
  if fixedh!=None: xx[2]=fixedh;bounds[2]=(fixedh,fixedh)
  res=minimize(NLL,xx*condition,args=(cases[place],vocnum[place],sig0,asc,precases[prereduce(place)]),bounds=bounds*np.repeat(condition,2).reshape([len(bounds),2]),method="SLSQP",options=minopts)
  if not res.success:
    print(res)
    print(place)
    print("xx =",xx)
    print("lcases =",list(cases[place]))
    print("lprecases =",precases[prereduce(place)])
    print("lvocnum =",vocnum[place])
    for x in sorted(list(opts)): print("%s:"%x,opts[x])
    print("bounds =",bounds)
    print("nweeks, ndays, minday, lastweek =",nweeks,",",ndays,",",minday,",",lastweek)
    raise RuntimeError(res.message)
  xx=res.x/condition

  # Work out the useful constant terms. They don't affect the optimisation, but they allow optimisation over hyperparameters
  
  
  # If 'statphase', make the log likelihood a better approximation to log(integral over all parameters) using stationary phase approximation
  if statphase:
    N=bmN+4
    eps=1e-3
    H=np.zeros([N,N])
    for i in range(N-1):
      for j in range(i+1,N):
        v=0
        eps1=eps/condition[i]
        eps2=eps/condition[j]
        for (s1,s2) in [(-1,-1),(-1,1),(1,-1),(1,1)]:
          x=np.copy(xx)
          x[i]+=s1*eps1
          x[j]+=s2*eps2
          v+=s1*s2*NLL(x*condition,cases[place],vocnum[place],sig0,asc,precases[prereduce(place)])
        e=v/(4*eps1*eps2)
        H[i,j]=e
        H[j,i]=e
    for i in range(N):
      x=np.copy(xx)
      v=0
      eps1=eps/condition[i]
      for s in [-1,0,1]:
        x=np.copy(xx)
        x[i]+=s*eps1
        v+=(s*s*3-2)*NLL(x*condition,cases[place],vocnum[place],sig0,asc,precases[prereduce(place)])
      H[i,i]=v/eps1**2
    det=np.linalg.det(H)
    if det<=0: print("Warning: Hessian not positive for %s. Can't make corrected log likelihood."%place);det=1
  else:
    det=1
    
  return res.x/condition,res.fun+log(det)/2

def printplaceinfo(place):
  print(place)
  print("="*len(place))
  print()
  print("                        Nonvar    Var   Seen")
  for w in range(nweeks):
    day0,day1=lastweek-(nweeks-w)*voclen+1,lastweek-(nweeks-1-w)*voclen
    print(daytodate(day0),"-",daytodate(day1),"%6d %6d %6.0f"%(vocnum[place][w][0],vocnum[place][w][1],sum(cases[place][day0-minday:day1-minday+1])))
  print()

def fullprint(AA,BB,lvocnum,lcases,h0,Tmin=None,Tmax=None,area=None):
  print("A      = estimated number of new cases of non-B.1.617.2 on this day multiplied by the ascertainment rate")
  print("B      = estimated number of new cases of B.1.617.2 on this day multiplied by the ascertainment rate")
  print("Pred   = predicted number of cases seen this day = A+B")
  print("Seen   = number of cases seen this day, after weekday adjustment")
  print("PredV1 = p*Pred, where p = proportion of non-B.1.617.2 amongst variant counts")
  print("PredV2 = (1-p)*Pred")
  print("SeenV1 = p*Seen")
  print("SeenV2 = (1-p)*Seen")
  print("Q      = estimated reproduction rate of non-B.1.617.2 on this day")
  print("R      = estimated reproduction rate of B.1.617.2 on this day")
  print()
  # Need the extra decimal places to make graphs look smooth
  if area!=None and args.graph_filename!=None:
    graphdata=args.graph_filename+'_'+area+'.dat'
    graphfp=open(graphdata,'w')
  else:
    graphfp=None
  def mprint(*a,**b):
    print(*a,**b)
    if graphfp!=None: print(*a,**b,file=graphfp)
  mprint("#     Date         A         B      Pred      Seen      PredV1    PredV2    SeenV1    SeenV2          Q       R")
  for i in range(ndays):
    day=minday+i
    pred,seen=asc*(AA[i]+BB[i]),lcases[i]
    mprint(daytodate(day),"%9.2f %9.2f %9.2f %9.2f"%(asc*AA[i],asc*BB[i],pred,seen),end='')
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks and lvocnum[week].sum()>0:
      p=lvocnum[week][0]/lvocnum[week].sum()
      mprint("   %9.2f %9.2f %9.2f %9.2f "%(p*pred,(1-p)*pred,p*seen,(1-p)*seen),end='')
    else:
      mprint("           -         -         -         - ",end='')
    if i<ndays-1:
      g=log(AA[i+1]/AA[i])
      Q,R=(exp(g*mgt),exp((g+h0)*mgt))
      mprint("   %7.4f %7.4f"%(Q,R))
    else:
      mprint()
  T=(R/Q-1)*100
  ETA="Estimated transmission advantage: %.0f%%"%T
  if Tmin!=None: ETA+=" (%.0f%% - %.0f%%), CI excludes model uncertainty"%(Tmin,Tmax)
  print(ETA)
  print()
  if graphfp!=None:
    graphfp.close()
    if area=="all":
      areacoveredhere=areacovered
    else:
      areacoveredhere=area
    # 'As of %s, estimated R(non-B.1.617.2)=%.2f R(B.1.617.2)=%.2f\\n'%(daytodate(minday+ndays-3),Q,R)+
    for yaxis in ["lin","log"]:
      graphfn=args.graph_filename+'_'+areacoveredhere+'_'+yaxis+'.png'
      po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
      # Use this write function to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
      def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))
      write('set xdata time')
      write('set key top left')
      write('set timefmt "%Y-%m-%d"')
      write('set format x "%Y-%m-%d"')
      write('set xtics nomirror rotate by 45 right offset 0.5,0')
      write('set label "As of %s:\\nEstimated R(non-B.1.617.2)=%.2f\\nEstimated R(B.1.617.2)=%.2f\\n%s" at screen 0.48,0.9'%(daytodate(minday+ndays-3),Q,R,ETA))
      write('set terminal pngcairo font "sans,13" size 1920,1280')
      write('set bmargin 7;set lmargin 13;set rmargin 13;set tmargin 5')
      write('set output "%s"'%graphfn)
      write('set ylabel "New cases per day (scaled down to match ascertainment rate of %0.f%%)"'%(100*asc))
      if yaxis=="log": write('set logscale y')
      write('set title "Estimated new cases per day of non-B.1.617.2 and B.1.617.2 in %s\\n'%areacoveredhere+
            'Fit made on 2021-05-29 using https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/vocfit.py\\n'+
            'Data sources: %s, Government coronavirus api/dashboard"'%fullsource)
      write('plot "%s" u 1:2 with lines lw 3 title "Modelled non-B.1.617.2", "%s" u 1:3 with lines lw 3 title "Modelled B.1.617.2", "%s" u 1:4 with lines lw 3 title "Modelled total", "%s" u 1:5 with lines lt 6 lw 3 title "Confirmed cases (all variants, weekday adjustment)", "%s" u 1:6 lt 1 pt 6 lw 3 title "Proportion of non-B.1.617.2 scaled up to modelled total", "%s" u 1:7 lt 2 pt 6 lw 3 title "Proportion of B.1.617.2 scaled up to modelled total"'%((graphdata,)*6))
      p.close();po.wait()
      print("Written graph to %s"%graphfn)
  return Q,R

def printsummary(summary):
  print("Location                       Q     R      T")
  for place in places:
    (Q,R,T,Tmin,Tmax)=summary[place]
    print("%-25s  %5.2f %5.2f  %4.0f%%"%(place,Q,R,T),end='')
    if Tmin!=None: print(" ( %4.0f%% - %4.0f%% )"%(Tmin,Tmax))
    else: print()
  print()
  print("Q = point estimate of reproduction rate of non-B.1.617.2 on",daytodate(maxday-1))
  print("R = point estimate of reproduction rate of B.1.617.2 on",daytodate(maxday-1))
  print("T = estimated transmission advantage = R/Q as a percentage increase")
  print()

if mode=="local growth rates":
  summary={}
  for place in places:
    printplaceinfo(place)
    xx0,L0=optimiseplace(place)
    if 0:
      print("Log likelihood",-L0)
      xx1,L1=optimiseplace(place,statphase=True)
      print("Corrected log likelihood",-L1)
    AA,BB,GG=expand(xx0)
    h0=xx0[2]
    ff=[0,L0,0]
    eps=0.01
    for i in [-1,1]:
      xx,L=optimiseplace(place,hint=xx0,fixedh=h0+i*eps)
      ff[i+1]=L
    # Use observed Fisher information to make confidence interval
    fi=(ff[0]-2*ff[1]+ff[2])/eps**2
    if fi>0:
      dh=1.96/sqrt(fi)
      (Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-dh,h0,h0+dh]]
    else:
      (Tmin,T,Tmax)=[None,(exp(h0*mgt)-1)*100,None]
    print("Locally optimised growth advantage")
    Q,R=fullprint(AA,BB,vocnum[place],cases[place],h0,Tmin,Tmax,area=place)
    summary[place]=(Q,R,T,Tmin,Tmax)
  print()
  printsummary(summary)

if type(mode)==tuple and mode[0]=="fixed growth rate":
  summary={}
  for place in places:
    printplaceinfo(place)
    h0=mode[1]
    xx0,L0=optimiseplace(place,fixedh=h0)
    AA,BB,GG=expand(xx0)
    T=(exp(h0*mgt)-1)*100
    print("Predetermined growth advantage")
    Q,R=fullprint(AA,BB,vocnum[place],cases[place],h0,area=place)
    summary[place]=(Q,R,T,None,None)
  print()
  printsummary(summary)
  
if mode=="global growth rate":
  ndiv=11;hmin=0.03;hmax=0.15
  logp=np.zeros(ndiv)
  for place in places:
    printplaceinfo(place)
    print("    h     T    log lik")
    xx,L0=optimiseplace(place)
    for i in range(ndiv):
      h=(hmin+(hmax-hmin)*i/(ndiv-1))
      xx,L=optimiseplace(place,hint=xx,fixedh=h)
      logp[i]+=L0-L
      print("%5.3f %5.3f  %9.2f"%(h,exp(h*mgt),logp[i]))
    print()
    sys.stdout.flush()
  print()
  i=np.argmax(logp)
  if i==0 or i==ndiv-1:
    print("Can't properly estimate best transmission factor or confidence interval because the maximum is at the end")
    imax=i
    c=0.1
  else:
    b=(logp[i+1]-logp[i-1])/2
    c=2*logp[i]-(logp[i+1]+logp[i-1])
    imax=i+b/c
  irange=1.96/sqrt(c)
  h0=(hmin+(hmax-hmin)*imax/(ndiv-1))
  dh=(hmax-hmin)*irange/(ndiv-1)
  (Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-dh,h0,h0+dh]]
  print("Combined growth advantage per day: %.3f (%.3f - %.3f)"%(h0,h0-dh,h0+dh))
  print("Combined transmission advantage: %.0f%% (%.0f%% - %.0f%%) (assuming fixed generation time of %g days)"%(T,Tmin,Tmax,mgt))
  print()
  print("Re-running using global optimum growth advantage",h0)
  print()

  summary={}
  TAA=np.zeros(ndays)
  TBB=np.zeros(ndays)
  for place in places:
    printplaceinfo(place)
    xx0,L0=optimiseplace(place,fixedh=h0)
    AA,BB,GG=expand(xx0)
    TAA+=AA;TBB+=BB
    T=(exp(h0*mgt)-1)*100
    print("Globally optimised growth advantage")
    Q,R=fullprint(AA,BB,vocnum[place],cases[place],h0)
    summary[place]=(Q,R,T,None,None)
  print()
  printsummary(summary)
  
  print("Total predicted counts using global optimum growth advantage")
  print()

  print("Combined results using globally optimised growth advantage")
  Q,R=fullprint(TAA,TBB,sum(vocnum.values()),[sum(cases[place][i] for place in places) for i in range(ndays)],h0,Tmin,Tmax,area="all")
  
  (Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-dh,h0,h0+dh]]
  print("Combined growth advantage per day: %.3f (%.3f - %.3f)"%(h0,h0-dh,h0+dh))
  print("Combined transmission advantage: %.0f%% (%.0f%% - %.0f%%) (assuming fixed generation time of %g days)"%(T,Tmin,Tmax,mgt))
