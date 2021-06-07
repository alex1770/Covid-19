from stuff import *
import sys,re,argparse,pickle
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

# LTLA_age.csv from https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesBySpecimenDateAgeDemographics&format=csv

def sanitise(fn): return fn.replace(' ','_').replace("'","")

ltlaengdata=loadcsv("Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")
ltla2region=dict(zip(ltlaengdata['LAD20CD'],ltlaengdata['RGN20NM']))
ltla2name=dict(zip(ltlaengdata['LAD20CD'],map(sanitise,ltlaengdata['LAD20NM'])))
ltlapopdata=loadcsv("LTLA-NIMS-populations.csv")

# 1. Sanger uses LAD19 except E06000053 (Isles of Scilly) isn't present. I assume it's fused into E06000052 (Cornwall).
# 2. Dashboard/api use LAD19 except that E09000001 (City of London) has been fused into E09000012 (Hackney).
#    and E06000053 (Isles of Scilly) is fused into E06000052 (Cornwall).
# 3. NIMS (population and vaccine data) uses LAD20 (E070000[4567] fused into E06000060).
# Easier to fuse than to split, so standardise on LAD20 with E09000001 fused into E09000012 and E06000053 fused into E06000052.
fuseltla=dict(zip(ltlaengdata['LAD20CD'],ltlaengdata['LAD20CD']))
fuseltla['E07000004']='E06000060'
fuseltla['E07000005']='E06000060'
fuseltla['E07000006']='E06000060'
fuseltla['E07000007']='E06000060'
fuseltla['E09000001']='E09000012'
fuseltla['E06000053']='E06000052'

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
### End Model ###


### Options ###

source="Sanger"
#source="COG-UK"
#source="SGTF"

# Number of parameters to optimise
#N=4;r0=0.3
N=3

vaxeffecttime=20# Days before vaccine is presumed to have a decent effect

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

if source=="Sanger":
  fullsource="Wellcome Sanger Institute"
  assert voclen==7
  sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
  
  lastweek=datetoday(max(sanger['WeekEndDate']))
  nweeks=(lastweek-firstweek)//voclen+1
  # Sanger week number is nweeks-1-(lastweek-day)//voclen

  # Get Sanger (variant) data into a suitable form
  vocnum={}#ltla: np.zeros([nweeks,2],dtype=int) for ltla in ltla2region if ltla[0]=='E'}
  rvocnum={}
  for (date,lad19,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
    day=datetoday(date)
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks:
      ltla=fuseltla[lad19]
      if ltla in ltlaexclude or not includeltla(ltla,ltlaset): continue
      if ltla not in vocnum: vocnum[ltla]=np.zeros([nweeks,2],dtype=int)
      if var=="B.1.617.2": vocnum[ltla][week][1]+=n
      else: vocnum[ltla][week][0]+=n
      place=ltla2region[ltla]
      if place not in rvocnum: rvocnum[place]=np.zeros([nweeks,2],dtype=int)
      if var=="B.1.617.2": rvocnum[place][week][1]+=n
      else: rvocnum[place][week][0]+=n
elif source=="COG-UK":
  fullsource="COG-UK"
  cog=loadcsv("cog_metadata.csv")
  lastweek=datetoday(max(cog['sample_date']))-discardcogdays
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

tvocnum=sum(rvocnum.values())

places=sorted([x for x in set(fuseltla.values()) if x not in ltlaexclude and includeltla(x,ltlaset)])

# Restrict to places for which there is at least some of each variant
#okplaces=set([place for place in places if place in vocnum and vocnum[place][:,0].sum()>0 and vocnum[place][:,1].sum()>0])
okplaces=set(places)

places=list(okplaces)
places.sort()# Alphabetical order

ltlapop={}
regionpop={}
for (lad20,n1,n2) in zip(ltlapopdata['LTLA Code'],ltlapopdata['Under 16'],ltlapopdata['16+']):
  ltla=fuseltla[lad20]
  ltlapop[ltla]=n1+n2
  region=ltla2region[ltla]
  regionpop[region]=regionpop.get(region,0)+n1+n2
rpopratio={ltla:ltlapop[ltla]/regionpop[ltla2region[ltla]] for ltla in ltlapop}
tpop=sum(regionpop.values())
tpopratio={ltla:ltlapop[ltla]/tpop for ltla in ltlapop}

# Convert daily growth rate & uncertainty into R-number-based description
# dh = 1 standard deviation
def Rdesc(h0,dh):
  (Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-zconf*dh,h0,h0+zconf*dh]]
  return "%.0f%% (%.0f%% - %.0f%%)"%(T,Tmin,Tmax)

# Need to scale the variables being optimised over to keep SLSQP happy
condition=np.zeros(N)+1

# Up to 5 parameters to be optimised:
# 0: g    Difference of weekly growth rates (g(V1)-g(V0)) for unvaccinated people
# 1: rw   Weight of region counts in hierarchy
# 2: tw   Weight of total counts in hierarchy
# 3: r1   RR (= 1-VE) of vaccine for variant 1 (delta)
# 4: r0   RR (= 1-VE) of vaccine for variant 0 (alpha)

# Return negative log likelihood (negative because scipy can only minimise, not maximise)
# If const is true then add in all the constant terms (that don't affect the optimisation)
def NLL(xx_conditioned,const=False,pic=False):
  xx=xx_conditioned/condition
  tot=0

  # Prior on G
  #a0=log(lcases[0]+.5)
  #tot+=-((xx[0]-a0)*isd0)**2/2
  #if const: tot-=log(2*pi/isd0**2)/2
  
  # Prior on rweight
  #tot+=-((xx[1]-(a0-4))*isd1)**2/2
  #if const: tot-=log(2*pi/isd1**2)/2
  
  # Prior on tweight
  #tot+=-((xx[2]-(a0-4))*isd1)**2/2
  #if const: tot-=log(2*pi/isd1**2)/2

  # Prior on r1
  #tot+=-(xx[3]*isd2)**2/2
  #if const: tot-=log(2*pi/isd2**2)/2
  
  # Prior on r0
  #tot+=-(xx[4]*isd2)**2/2
  #if const: tot-=log(2*pi/isd2**2)/2

  if pic: fp=open("temp","w")
  for place in places:
    if place not in vocnum: continue
    vv=vocnum[place]
    rv=rvocnum[ltla2region[place]]
    rpr=rpopratio[place]
    tpr=tpopratio[place]
    for w in range(nweeks-1):
      rho=exp(xx[0])
      AB=vv[w]+xx[1]*rpr*rv[w]+xx[2]*tpr*tvocnum[w]
      CD=vv[w+1]
      if N>=4:
        p=pvax[place][w]
        r1=xx[3]
        if N>=5: r0_=xx[4]
        else: r0_=r0
        rho*=(1-p+p*r1)/(1-p+p*r0_)
        if pic: print("%s.%d  %8.5f  %12g  %12g"%(place,w,p,CD[0]/AB[0],CD[1]/AB[1]),file=fp)
      s=AB[0]+rho*AB[1]
      tot+=CD[0]*log(AB[0]/s)+CD[1]*log(rho*AB[1]/s)
      if const: tot+=gammaln(CD[0]+CD[1]+1)-gammaln(CD[0]+1)-gammaln(CD[1]+1)# Could make a table
  if pic: fp.close()
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
def optimise(hint=[0.7,.25,.2,0.3,0.3][:N],statphase=False):
  xx=np.copy(hint)
  bounds=[(-5,10),(1e-2,100),(1e-2,100),(0.01,5),(0.01,5)][:N]
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

vaxdir='VaccinationData'
vaxdat={}
for x in sorted(os.listdir(vaxdir)):
  f=x.find('20')
  vaxdat[x[f:f+10]]=loadcsv(os.path.join(vaxdir,x))

def parseage_api(age):
  if '_' in age: x=age.split('_');a0,a1=int(x[0]),int(x[1])+1
  elif age[-1]=='+': a0,a1=int(age[:-1]),150
  elif age=='unassigned': a0,a1=0,150
  else: raise RuntimeError("Unrecognised age band: "+age)
  return (a0,a1)
  
def parseage_nims(age):
  return (a0,a1)
  
ltlaagecachedir="LTLA_age_cache"
ltlaagecachename="LTLA_age_weekly_%s_+%dweeks.pickle"%(daytodate(firstweek),nweeks)
fn=os.path.join(ltlaagecachedir,ltlaagecachename)
os.makedirs(ltlaagecachedir,exist_ok=True)
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    caseages=pickle.load(fp)
else:
  la=loadcsv("LTLA_age.csv")
  date0=daytodate(lastweek)
  date1=max(la['date'])
  if date1<date0: raise RuntimeError("LTLA_age.csv not up to date. Need up to %s but it ends at %s."%(date0,date1))
  caseages={}
  for lad19,date,age,cases in zip(la['areaCode'],la['date'],la['age'],la['cases']):
    ltla=fuseltla[lad19]
    day=datetoday(date)
    w=nweeks-1-(lastweek-day)//voclen
    if w>=0 and w<nweeks:
      if ltla not in caseages: caseages[ltla]={}
      a=parseage_api(age)
      if a!=(0,60) and a!=(60,150):
        caseages[ltla].setdefault(a,[0]*nweeks)[w]+=cases
  with open(fn,'wb') as fp:
    pickle.dump(caseages,fp)

    
from random import random,seed
seed(42)#alter
pvax={place:[0]*(nweeks-1) for place in places}
for w in range(nweeks-1):
  date=daytodate(firstweek+w*7+10-vaxeffecttime)
  id=max(dt for dt in vaxdat if dt<date)
  print("Using NIMS vax w/e %s to correspond to Sanger w/e %s -> w/e %s"%(id,daytodate(firstweek+w*7),daytodate(firstweek+w*7+7)))
  v=vaxdat[id]
  for (i,place) in enumerate(v['LTLA Code']):
    if place in places:
      n=sum(v[key][i] for key in v if key[0]=='D')
      pvax[place][w]=min(n/ltlapop[place],1)

xx,L=optimise()
print("Variables:",xx)
print("Log likelihood:",L)
NLL(xx*condition,const=True,pic=True)
H=Hessian(xx)
h=xx[0]/7;dh=1/sqrt(H[0,0])/7
print("Logarithmic growth rate advantage/day: %.1f%% (%.1f%% - %.1f%%)"%(h*100,(h-zconf*dh)*100,(h+zconf*dh)*100))
print("Multiplicative growth rate advantage/day: %.1f%% (%.1f%% - %.1f%%)"%((exp(h)-1)*100,(exp(h-zconf*dh)-1)*100,(exp(h+zconf*dh)-1)*100))
print("R-number advantage: %.2f (%.2f - %.2f)"%(exp(mgt*h),exp(mgt*(h-zconf*dh)),exp(mgt*(h+zconf*dh))))