from stuff import *
import sys,re,argparse
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.special import gammaln
from math import log,exp,sqrt,sin,pi
import numpy as np
from subprocess import Popen,PIPE
from datetime import datetime

# (Make it auto download files?)
# Get ltla.csv from https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv
# Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv
# COG-UK data from https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv
# SGTF   data from Fig.16 Tech Briefing 12: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/988608/Variants_of_Concern_Technical_Briefing_12_Data_England.xlsx

def sanitise(fn): return fn.replace(' ','_').replace("'","")

apicases=loadcsv("ltla.csv")
ltlaengdata=loadcsv("Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv")
ltlaukdata=loadcsv("Local_Authority_District_to_Country_(April_2019)_Lookup_in_the_United_Kingdom.csv")
ltla2ltla=dict(zip(ltlaukdata['LAD19CD'],ltlaukdata['LAD19CD']))
ltla2uk=dict((ltla,"UK") for ltla in ltlaukdata['LAD19CD'])
ltla2country=dict(zip(ltlaukdata['LAD19CD'],ltlaukdata['CTRY19NM']))
ltla2region=dict(ltla2country,**dict(zip(ltlaengdata['LAD19CD'],ltlaengdata['RGN19NM'])))
ltla2name=dict(zip(ltlaukdata['LAD19CD'],map(sanitise,ltlaukdata['LAD19NM'])))
ltla2pop=loadcsv("LTLA-NIMS-populations.csv")

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
  elif ltlaset=="test":
    return ltla2region[ltla]=='London' and ltla<'E09000010'
  elif ltlaset=="Bolton":
    return ltla=='E08000001'
  elif ltlaset=="Hartlepool":
    return ltla=='E06000001'
  elif ltlaset=="NE":
    return ltla2region[ltla]=='North East'
  elif ltlaset=="All":
    return True
  else:
    raise RuntimeError("Unrecognised ltla set "+ltlaset)

parser=argparse.ArgumentParser()
parser.add_argument('-l', '--load-options',   help='Load options from a file')
parser.add_argument('-s', '--save-options',   help='Save options to a file')
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
# n_i ~ NB(mean=p(A_i+B_i),var=mean/nif1)
# r_j ~ BetaBinomial(r_j+s_j, A_{I_j}nif2/(1-nif2), B_{I_j}nif2/(1-nif2))  (A_{I_j} means sum_{i in I_j}A_i)
# g_0 ~ N(g_{-1},v_{-1})
# X_n ~ N(0,1)
# N=ndays-1
# g_i = g_0 + bmscale*sqrt(N)*(i/N*X_0 + sqrt(2)/pi*sum_n e^{-(n*bmsig/N)^2/2}sin(n*pi*i/N)*X_n/n)
# h ~ N(0,tau^2)
#
### End Model ###


### Options ###

source="Sanger"
#source="COG-UK"
#source="SGTF"

N=2# Number of parameters to optimise

# Can choose location size from "LTLA", "region", "country", "UK"
# Sanger works with LTLA, region, country
# COG-UK works with country, UK
# SGTF works with region, country
locationsize="LTLA"

ltlaexclude=set()
# ltlaexclude=set(['E08000001'])# This would exclude Bolton
ltlaset="All"
#ltlaset="London"
#ltlaset="Bolton"
#ltlaset="Hartlepool"

# Will plot graph of these locations even if only encountered during subdivision of global growth mode
specialinterest=set(['E08000001'])

mgt=5# Mean generation time in days

# Earliest day to use case data
minday=datetoday('2021-04-01')# Inclusive

# Earliest day to use VOC count data, given as end-of-week. Will be rounded up to match same day of week as lastweek.
#firstweek=minday+6
firstweek=datetoday('2021-05-01')

nif1=0.048 # Non-independence factor (1/overdispersion) for cases (less than 1 means information is downweighted)
nif2=0.255 # Non-independence factor (1/overdispersion) for VOC counts (ditto)
isd0=1.0   # Inverse sd for prior on starting number of cases of non-B.1.617.2: assume starts off similar to total number of cases
isd1=0.3   # Inverse sd for prior on starting number of cases of B.1.617.2 (0.3 is very weak)
isd2=1     # Inverse sd for prior on transmission advantage (as growth rate per day). 0 means uniform prior. 1 is very weak.

# Prior linking initial daily growth rate to estimate from pre-B.1.617.2 era
sig0=0.004

# Timescale in days over which growth rate can change significantly
# (lower = more wiggles)
bmsig=25

# Lengthscale for filtered Brownian motion
# (higher = greater amplitude for the wiggles)
bmscale=0.01

# Case ascertainment rate
asc=0.4

# Discard this many cases at the end of the list of cases by specimen day
discardcasedays=3# Pro tem to allow for Wales and NI late reporting at weekends and bank holidays

# Discard this many days of the latest COG data
discardcogdays=2

minopts={"maxiter":10000,"eps":1e-5}

mode="local growth rates"
#mode="global growth rate"
#mode="fixed growth rate",0.1

voclen=(1 if source=="COG-UK" else 7)

conf=0.95
nsamp=2000

model="scaledpoisson"
#model="NBBB"
#model="NBBB+magicprior"

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
  "Inverse sd for prior on initial non-B.1.617.2": isd0,
  "Inverse sd for prior on initial B.1.617.2": isd1,
  "Inverse sd for prior on growth": isd2,
  "Sigma (prior on daily growth rate change)": sig0,
  "Timescale of growth rate change (days)": bmsig,
  "Lengthscale for filtered Brownian motion": bmscale,
  "Case ascertainment rate": asc,
  "Number of days of case data to discard": discardcasedays,
  "Number of days of COG-UK data to discard": discardcogdays,
  "Minimiser options": minopts,
  "Length of time period over which VOC counts are given (days)": voclen,
  "Confidence level": conf,
  "Number of samples for confidence calcultions in hierachical mode": nsamp,
  "Model": model
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
isd0=opts["Inverse sd for prior on initial non-B.1.617.2"]
isd1=opts["Inverse sd for prior on initial B.1.617.2"]
isd2=opts["Inverse sd for prior on growth"]
sig0=opts["Sigma (prior on daily growth rate change)"]
bmsig=opts["Timescale of growth rate change (days)"]
bmscale=opts["Lengthscale for filtered Brownian motion"]
asc=opts["Case ascertainment rate"]
discardcasedays=opts["Number of days of case data to discard"]
discardcogdays=opts["Number of days of COG-UK data to discard"]
minopts=opts["Minimiser options"]
voclen=opts["Length of time period over which VOC counts are given (days)"]
conf=opts["Confidence level"]
nsamp=opts["Number of samples for confidence calcultions in hierachical mode"]
model=opts["Model"]

zconf=norm.ppf((1+conf)/2)

print("Options:")
print()
for x in sorted(list(opts)): print("%s:"%x,opts[x])
print()
sys.stdout.flush()

np.set_printoptions(precision=3,linewidth=150)

if ltlaset=="All":
  if source=="COG-UK": areacovered="UK"
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
  vocnum={ltla: np.zeros([nweeks,2],dtype=int) for ltla in ltla2region if ltla[0]=='E'}
  rvocnum={}
  for (date,ltla,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
    if ltla in ltlaexclude or not includeltla(ltla,ltlaset): continue
    day=datetoday(date)
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks:
      place=ltla
      if var=="B.1.617.2": vocnum[place][week][1]+=n
      else: vocnum[place][week][0]+=n
      place=ltla2region[ltla]
      if place not in rvocnum: rvocnum[place]=np.zeros([nweeks,2],dtype=int)
      if var=="B.1.617.2": rvocnum[place][week][1]+=n
      else: rvocnum[place][week][0]+=n
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

places=sorted(list(set(sanger['LTLA'])))

# Restrict to places for which there is at least some of each variant
#okplaces=set([place for place in places if vocnum[place][:,0].sum()>0 and vocnum[place][:,1].sum()>0])
okplaces=set(places)

places=list(okplaces)
places.sort()# Alphabetical order

# Adjustment because Sanger and the api use LAD19 and NIMS uses LAD20 (or something)
# Old LTLA  Population      New LTLA  Population
# E07000004     215000     E06000060      581000
# E07000005     103000
# E07000006      75000
# E07000007     188000
ltlapop={}
regionpop={}
ltlapop['E07000004']=215000
ltlapop['E07000005']=103000
ltlapop['E07000006']= 75000
ltlapop['E07000007']=188000
regionpop[ltla2region['E07000004']]=581000
for (ltla,n1,n2) in zip(ltla2pop['LTLA Code'],ltla2pop['Under 16'],ltla2pop['16+']):
  if ltla!='E06000060':
    ltlapop[ltla]=n1+n2
    region=ltla2region[ltla]
    regionpop[region]=regionpop.get(region,0)+n1+n2
popratio={ltla:ltlapop[ltla]/regionpop[ltla2region[ltla]] for ltla in ltlapop}

# Convert daily growth rate & uncertainty into R-number-based description
# dh = 1 standard deviation
def Rdesc(h0,dh):
  (Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-zconf*dh,h0,h0+zconf*dh]]
  return "%.0f%% (%.0f%% - %.0f%%)"%(T,Tmin,Tmax)

# L=ndays-1
# i is time from start, in days
# t=i/L
# growth[i] = bmscale*sqrt(L)*(t*X_0 + sum_{n=1}^{N-1} sqrt(2)/pi*exp(-(n*bmsig/L)^2/2)*sin(n*pi*t)/n*X_n)
# where X_n ~ N(0,1),  n=0,...,N; N=ceil(4*L/bmsig), say

bmL=ndays-1
bmN=int(3*bmL/bmsig+1)
bmsin=[sin(r*pi/bmL) for r in range(2*bmL)]
bmweight=[0]+[sqrt(2)/pi*exp(-(n*bmsig/bmL)**2/2)/n for n in range(1,bmN)]
bmsin2=[np.array([bmweight[n]*bmsin[(i*n)%(2*bmL)] for n in range(bmN)]) for i in range(ndays)]

# Need to scale the variables being optimised over to keep SLSQP happy
condition=np.zeros(N)+1

# Up to 4 parameters to be optimised:
# 0: g    Difference of weekly growth rates (g(V1)-g(V0)) for unvaccinated people
# 1: w    Weight of region counts in hierarchy
# 2: r1   RR (= 1-VE) of vaccine for variant 1 (delta)
# 3: r0   RR (= 1-VE) of vaccine for variant 0 (alpha)

# Return negative log likelihood (negative because scipy can only minimise, not maximise)
# If const is true then add in all the constant terms (that don't affect the optimisation)
def NLL(xx_conditioned,const=False):
  xx=xx_conditioned/condition
  tot=0

  # Prior on G
  #a0=log(lcases[0]+.5)
  #tot+=-((xx[0]-a0)*isd0)**2/2
  #if const: tot-=log(2*pi/isd0**2)/2
  
  # Prior on w
  #tot+=-((xx[1]-(a0-4))*isd1)**2/2
  #if const: tot-=log(2*pi/isd1**2)/2
  
  # Prior on r1
  #tot+=-(xx[2]*isd2)**2/2
  #if const: tot-=log(2*pi/isd2**2)/2
  
  # Prior on r0
  #tot+=-(xx[2]*isd2)**2/2
  #if const: tot-=log(2*pi/isd2**2)/2

  for place in places:
    vv=vocnum[place]
    rv=rvocnum[ltla2region[place]]
    pr=popratio[place]
    for w in range(nweeks-1):
      rho=exp(xx[0])
      AB=vv[w]+xx[1]*pr*rv[w]
      CD=vv[w+1]
      s=AB[0]+rho*AB[1]
      tot+=CD[0]*log(AB[0]/s)+CD[1]*log(rho*AB[1]/s)
      if const: tot+=gammaln(CD[0]+CD[1]+1)-gammaln(CD[0]+1)-gammaln(CD[1]+1)# Could make a table

  return -tot

def Hessian(xx):
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
        v+=s1*s2*NLL(x*condition)
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
      v+=(s*s*3-2)*NLL(x*condition)
    H[i,i]=v/eps1**2
  return H

# Returns log likelihood
def optimise(hint=[0.,1,1,1][:N],statphase=False):
  xx=np.copy(hint)
  # bounds[2][0]=0 prejudges B.1.617.2 as being at least as transmissible as B.1.1.7. This helps SLSQP not get stuck in some cases
  # though would need to relax this constraint if dealing with other variants where it might not be true.
  bounds=[(-5,10),(1e-2,100),(0.01,5),(0.01,5)][:N]
  res=minimize(NLL,xx*condition,bounds=bounds,method="SLSQP",options=minopts)
  if not res.success:
    print(res)
    print(place)
    print("xx =",xx)
    for x in sorted(list(opts)): print("%s:"%x,opts[x])
    print("bounds =",bounds)
    raise RuntimeError(res.message)
  xx=res.x/condition

  # Work out log likelihood including constant terms
  LL=-NLL(res.x,const=True)
  
  # If 'statphase', make the log likelihood a better approximation to log(integral over all parameters) using stationary phase approximation
  if statphase:
    H=Hessian(xx)
    det=np.linalg.det(H)
    if det<=0: print("Warning: Hessian not positive for %s. Can't make corrected log likelihood."%place);det=1
    LL+=N*log(2*pi)/2-log(det)/2

  # Return optimum xx log likelihood
  return res.x/condition,LL

def evalconfidence(place,xx0):
  H=Hessian(xx0,cases[place],vocnum[place],sig0,asc,precases[prereduce(place)])
  Hcond=H/condition/condition[:,None]
  eig=np.linalg.eigh(Hcond)
  # np.diag(np.matmul(np.matmul(np.transpose(eig[1]),Hcond),eig[1])) ~= eig[0]
  if not (eig[0]>0).all(): print("Hessian not +ve definite so can't do full confidence calculation");return None,None,None,None
  nsamp=10000
  N=bmN+4
  t=norm.rvs(size=[nsamp,N])# nsamp x N
  sd=eig[0]**(-.5)# N
  u=t*sd# nsamp x N
  samp_cond=np.matmul(u,np.transpose(eig[1]))# nsamp x N
  samp=samp_cond/condition
  QQ=[];RR=[];TT=[]
  for i in range(nsamp):
    dx=samp[i]
    xx=xx0+dx
    AA,BB,GG=expand(xx)
    QQ.append((AA[-1]/AA[-2])**mgt)
    RR.append((BB[-1]/BB[-2])**mgt)
    TT.append(exp(mgt*xx[2]))
  QQ.sort();RR.sort();TT.sort()
  n0=int((1-conf)/2*nsamp)
  n1=int((1+conf)/2*nsamp)
  print("R_{B.1.1.7}   %6.3f - %6.3f"%(QQ[n0],QQ[n1]))
  print("R_{B.1.617.2} %6.3f - %6.3f"%(RR[n0],RR[n1]))
  print("T             %5.1f%% - %5.1f%%"%((TT[n0]-1)*100,(TT[n1]-1)*100))
  return QQ[n0],QQ[n1],RR[n0],RR[n1]


def printplaceinfo(place,using=''):
  name=ltla2name.get(place,place)+using
  print()
  print(name)
  print("="*len(name))
  print()
  print("                        Nonvar    Var   Seen")
  for w in range(nweeks):
    day0,day1=lastweek-(nweeks-w)*voclen+1,lastweek-(nweeks-1-w)*voclen
    print(daytodate(day0),"-",daytodate(day1),"%6d %6d %6.0f"%(vocnum[place][w][0],vocnum[place][w][1],sum(cases[place][day0-minday:day1-minday+1])))
  print()

def fullprint(AA,BB,lvocnum,lcases,T,Tmin=None,Tmax=None,Qmin=None,Qmax=None,Rmin=None,Rmax=None,area=None,using=''):
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
  if area!=None and args.graph_filename!=None:
    graphdata=sanitise(args.graph_filename+'_'+area+'.dat')
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
      Q,R=((AA[i+1]/AA[i])**mgt,(BB[i+1]/BB[i])**mgt)
      mprint("   %7.4f %7.4f"%(Q,R))
    else:
      mprint()
  # Note that T is not 100(R/Q-1) here because AA, BB are derived from a sum of locations each of which has extra transm T,
  # but because of Simpson's paradox, that doesn't been the cross ratio of AAs and BBs is also T.
  EQ="Estimated R(non-B.1.617.2) = %.2f"%Q
  if Qmin!=None: EQ+=" (%.2f - %.2f)"%(Qmin,Qmax)
  ER="Estimated R(B.1.617.2)       = %.2f"%R
  if Rmin!=None: ER+=" (%.2f - %.2f)"%(Rmin,Rmax)
  ETA="Estimated transmission advantage = %.0f%%"%T
  if Tmin!=None: ETA+=" (%.0f%% - %.0f%%)"%(Tmin,Tmax)
  if area!=None: print("Summary");print(area+using)
  print(EQ)
  print(ER)
  print(ETA)
  if Tmin!=None: ETA+="\\n(CI shows within-model statistical uncertainty, not model uncertainty)"
  print()
  if graphfp!=None:
    graphfp.close()
    now=datetime.utcnow().strftime('%Y-%m-%d')
    for yaxis in ["lin","log"]:
      graphfn=sanitise(args.graph_filename+'_'+area+'_'+yaxis+'.png')
      po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
      # Use this write function to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
      def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))
      write('set xdata time')
      write('set key top left')
      write('set timefmt "%Y-%m-%d"')
      write('set format x "%Y-%m-%d"')
      write('set xtics nomirror rotate by 45 right offset 0.5,0')
      write('set label "Location: %s\\nAs of %s:\\n%s\\n%s\\n%s" at screen 0.48,0.9'%(area+using,daytodate(minday+ndays-2),EQ,ER,ETA))
      write('set terminal pngcairo font "sans,13" size 1920,1280')
      write('set bmargin 7;set lmargin 13;set rmargin 13;set tmargin 5')
      write('set output "%s"'%graphfn)
      write('set ylabel "New cases per day (scaled down to match ascertainment rate of %0.f%%)"'%(100*asc))
      if yaxis=="log": write('set logscale y')
      write('set title "Estimated new cases per day of non-B.1.617.2 and B.1.617.2 in %s\\n'%(area+using)+
            'Fit made on %s using https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/vocfit.py\\n'%now+
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

xx,L=optimise()
print(xx)
print("Log likelihood:",L)
H=Hessian(xx)
h=xx[0]/7;dh=1/sqrt(H[0,0])/7
print("Growth rate advantage/day: %.1f%% (%.1f%% - %.1f%%)"%(h*100,(h-zconf*dh)*100,(h+zconf*dh)*100))
print("R-number advantage %.2f (%.2f - %.2f):"%(exp(mgt*h),exp(mgt*(h-zconf*dh)),exp(mgt*(h+zconf*dh))))
