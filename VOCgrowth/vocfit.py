from stuff import *
import sys,re,argparse,pickle,pytz
from scipy.optimize import minimize
from scipy.stats import norm,binom,bernoulli
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

ltlaengdata=loadcsv("Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv")
ltlaukdata=loadcsv("Local_Authority_District_to_Country_(April_2019)_Lookup_in_the_United_Kingdom.csv")
ltla2ltla=dict(zip(ltlaukdata['LAD19CD'],ltlaukdata['LAD19CD']))
ltla2uk=dict((ltla,"UK") for ltla in ltlaukdata['LAD19CD'])
ltla2country=dict(zip(ltlaukdata['LAD19CD'],ltlaukdata['CTRY19NM']))
ltla2region=dict(ltla2country,**dict(zip(ltlaengdata['LAD19CD'],ltlaengdata['RGN19NM'])))
ltla2name=dict(zip(ltlaukdata['LAD19CD'],map(sanitise,ltlaukdata['LAD19NM'])))

#variantset=["B.1.617.2", "AY."];nonvariantset=["B.1.1.7"];variant="Delta";nonvariant="Alpha"
#variantset=["AY.4.2"];nonvariantset=[""];variant="AY.4.2";nonvariant="non-AY.4.2"
variantset=["BA.1", "BA.2"];nonvariantset=[""];variant="Omicron";nonvariant="non-Omicron"
relative=variant+'/'+nonvariant

# Set bounds for relative daily growth rate
(hmin,hmax)=(0.0,0.2)
if variant=="B.1.617.2": (hmin,hmax)=(0.03,0.15)
if variant=="AY.4.2": (hmin,hmax)=(0.01,0.04)
if variant=="Omicron": (hmin,hmax)=(0.3,0.5)

def varmatch(var,pattern):
  if pattern=="": return True
  if pattern[-1:]!='.': return var==pattern
  return var[:len(pattern)]==pattern

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
# r_j, s_j = Variant counts of non-Delta, Delta in j^th week
# I_j      = set of days (week) corresponding to VOC counts r_j, s_j
# Assume chance of sequencing a case is a totally free parameter, and optimise over it
#
# Unknown:
# h   = daily growth advantage of Delta over other variants
# X_n = Fourier coefficients controlling the growth of non-Delta
# A_0 = initial count of non-Delta
# B_0 = initial count of Delta
#
# Likelihood:
# A_{i+1}=e^{g_i}A_i
# B_{i+1}=e^{g_i+h}B_i
# n_i ~ NB(mean=p(A_i+B_i),var=mean/nif1)
# r_j ~ BetaBinomial(r_j+s_j, A_{I_j}nif2/(1-nif2), B_{I_j}nif2/(1-nif2))  (A_{I_j} means sum_{i in I_j}A_i)
# g_0 ~ N(g_{-1},v_{-1})
# X_n ~ N(0,1)
# L=ndays+2*bmsig
# g_i = g_0 + bmscale*sqrt(L)*(i/L*X_0 + sqrt(2)/pi*sum_n e^{-(n*bmsig/L)^2/2}sin(n*pi*i/L)*X_n/n)
# h ~ N(0,tau^2)
#
### End Model ###


### Options ###

#source="Sanger"
#source="COG-UK"
source="SGTF"

# Can choose location size from "LTLA", "region", "country", "UK"
# Sanger works with LTLA, region, country
# COG-UK works with country, UK, and partially with coglab (coglab doesn't match up with case-count areas, so can only do VOC-only growth rates)
# SGTF works with region, country
#locationsize="coglab"
#locationsize="LTLA"
locationsize="region"
#locationsize="country"
#locationsize="UK"

ltlaexclude=set()
#ltlaexclude=set(['E08000001','E12000002'])# Bolton, Manchester
#ltlaexclude=set(['E08000001','E12000002']+[x for x in ltla2region if ltla2region[x]=='London'])# Bolton, Manchester, London
ltlaset="All"
#ltlaset="London"
#ltlaset="Bolton"
#ltlaset="Hartlepool"

# Will plot graph of these locations even if only encountered during subdivision of global growth mode
specialinterest=set()#['E08000001'])

mgt=5# Mean generation time in days

# Earliest day to use case data
minday=datetoday('2021-11-10')# Inclusive

# Earliest day to use VOC count data, given as end-of-week. Will be rounded up to match same day of week as lastweek.
#firstweek=minday+6
firstweek=datetoday('2021-11-10')

nif1=0.048 # Non-independence factor (1/overdispersion) for cases (less than 1 means information is downweighted)
nif2=0.255 # Non-independence factor (1/overdispersion) for VOC counts (ditto)
isd0=1.0   # Inverse sd for prior on starting number of cases of non-Delta: assume starts off similar to total number of cases
isd1=0.3   # Inverse sd for prior on starting number of cases of Delta (0.3 is very weak)
isd2=1     # Inverse sd for prior on competitive advantage (as growth rate per day). 0 means uniform prior. 1 is very weak.

# Prior linking initial daily growth rate to estimate from pre-Delta era
sig0=0.004

# Timescale in days over which growth rate can change significantly
# (lower = more wiggles)
bmsig=25

# Lengthscale for filtered Brownian motion
# (higher = greater amplitude for the wiggles)
# bmscale will be set below - no longer user specifiable
bmscale=0.01

# Case ascertainment rate
asc=0.4

# Discard this many cases at the end of the list of cases by specimen day (may be increased later if Wales is in the mix)
discardcasedays=2

# Discard this many days of the latest COG data
discardcogdays=2

# Collect together all locations without positive entries into one combined "Other" location
# (Makes little difference in practice)
bundleremainder=True

minopts={"maxiter":10000,"eps":1e-4,'ftol':1e-12}

#mode="local growth rates"
mode="global growth rate"
#mode="fixed growth rate",0.1

voclen=(1 if source=="COG-UK" or source=="SGTF" else 7)

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
  "Inverse sd for prior on initial non-Delta": isd0,
  "Inverse sd for prior on initial Delta": isd1,
  "Inverse sd for prior on growth": isd2,
  "Sigma (prior on daily growth rate change)": sig0,
  "Timescale of growth rate change (days)": bmsig,
  "Lengthscale for filtered Brownian motion": bmscale,
  "Case ascertainment rate": asc,
  "Number of days of case data to discard": discardcasedays,
  "Number of days of COG-UK data to discard": discardcogdays,
  "Bundle remainder": bundleremainder,
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
isd0=opts["Inverse sd for prior on initial non-Delta"]
isd1=opts["Inverse sd for prior on initial Delta"]
isd2=opts["Inverse sd for prior on growth"]
sig0=opts["Sigma (prior on daily growth rate change)"]
bmsig=opts["Timescale of growth rate change (days)"]
bmscale=opts["Lengthscale for filtered Brownian motion"]
asc=opts["Case ascertainment rate"]
discardcasedays=opts["Number of days of case data to discard"]
discardcogdays=opts["Number of days of COG-UK data to discard"]
bundleremainder=opts["Bundle remainder"]
minopts=opts["Minimiser options"]
voclen=opts["Length of time period over which VOC counts are given (days)"]
conf=opts["Confidence level"]
nsamp=opts["Number of samples for confidence calcultions in hierachical mode"]
model=opts["Model"]

# Make sure this is at least 2 if Wales is in the mix, because it has later reporting
discardcasedays=max(int(locationsize=="UK" or (source=="COG-UK" and locationsize=="country"))*2,discardcasedays)
opts["Number of days of case data to discard"]=discardcasedays

# if locationsize=="LTLA": bmscale=0.1
# elif locationsize=="region": bmscale=0.03
# else: bmscale=0.01
# opts["Lengthscale for filtered Brownian motion"]=bmscale

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

apicases=loadcsv("ltla.csv")
d=datetime.now(pytz.timezone("Europe/London"))
today=Date(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<16+10/60: today-=1# Dashboard/api updates at 4pm UK time
if max(apicases['date'])<today-1:
  import requests
  url='https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv'
  response=requests.get(url, timeout=10)
  if not response.ok: raise RuntimeError(response.text)
  with open('ltla.csv','w') as fp: fp.write(response.text)
  apicases=loadcsv('ltla.csv')

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
      if any(varmatch(var,pat) for pat in variantset): vocnum[place][week][1]+=n
      elif any(varmatch(var,pat) for pat in nonvariantset): vocnum[place][week][0]+=n
elif source=="COG-UK":
  fullsource="COG-UK"
  cog=loadcsv("cog_metadata.csv")
  lastweek=datetoday(max(cog['sample_date']))-discardcogdays
  assert maxday>=lastweek
  nweeks=(lastweek-firstweek)//voclen+1
  # Week number is nweeks-1-(lastweek-day)//voclen
  
  if locationsize=="coglab":
    reduceltla=None;bundleremainder=False
    reducecog=coglab2coglab
  elif locationsize=="country":
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
      if any(varmatch(var,pat) for pat in variantset): vocnum[place][week][1]+=1
      elif any(varmatch(var,pat) for pat in nonvariantset): vocnum[place][week][0]+=1
elif source=="SGTF":
  assert voclen==1
  l=[x for x in os.listdir('.') if x[:19]=='sgtf_regionepicurve']
  if l==[]: raise RuntimeError("No sgtf_regionepicurve csv file found in current directory")
  sgtf=loadcsv(max(l))
  fullsource="SGTF data from Omicron daily overview, last specimen date "+max(sgtf['specimen_date'])
  lastweek=max(datetoday(x) for x in sgtf['specimen_date'])
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
  background=[0,0]
  for (date,region,var,n) in zip(sgtf['specimen_date'],sgtf['PHEC_name'],sgtf['sgtf'],sgtf['n']):
    day=datetoday(date)
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks:
      place=reducesgtf(region)
      if place not in vocnum: vocnum[place]=np.zeros([nweeks,2],dtype=int)
      vocnum[place][week][int("SGTF" in var)]+=n
    if date>=Date('2021-10-01') and date<Date('2021-11-10'):
      background[int("SGTF" in var)]+=n
  # Adjust for non-Omicron SGTFs, based on the assumption that these are in a non-location-dependent proportion to the number of non-Omicron cases
  f=background[1]/background[0]
  for place in vocnum:
    for week in range(nweeks):
      vocnum[place][week][1]=max(vocnum[place][week][1]-int(f*vocnum[place][week][0]+.5),0)
else:
  raise RuntimeError("Unrecognised source: "+source)

# speccasesadjust[d][r] = chance that a specimen from day of the week d (Monday=0) is reported (by) r days later
speccasesadjust=np.array([[ 0.    , 0.3954, 0.8823, 0.9664, 0.9914, 0.9975, 0.9989, 1.    ],
                          [ 0.    , 0.3264, 0.8731, 0.9707, 0.9893, 0.9962, 1.0003, 1.    ],
                          [ 0.    , 0.3561, 0.8745, 0.9604, 0.983 , 0.989 , 0.9936, 1.    ],
                          [ 0.    , 0.3511, 0.8699, 0.9567, 0.9768, 0.9873, 0.9974, 1.    ],
                          [ 0.    , 0.3132, 0.8363, 0.942 , 0.9701, 0.988 , 0.9964, 1.    ],
                          [ 0.    , 0.2943, 0.8404, 0.9218, 0.9595, 0.984 , 0.994 , 1.    ],
                          [ 0.    , 0.4585, 0.9093, 0.961 , 0.9859, 0.9934, 0.9987, 1.    ]])
infinity=speccasesadjust.shape[1]# Cases are considered stable at infinity-1 days after specimen taken
publishedday=datetoday(max(apicases['date']))+1# Day when api published results
monday=datetoday('2021-06-07')# Example of a Monday
specadj=np.zeros(ndays)
for d in range(ndays):
  day=minday+d
  r=publishedday-day;assert r>=1 and r<1000
  if r<infinity: specadj[d]=speccasesadjust[(day-monday)%7][r]
  else: specadj[d]=1
nif1a=nif1*specadj
lognif1a=np.log(nif1a)
log1mnif1a=np.log(1-nif1a)

# Simple weekday adjustment by dividing by the average count for that day of the week.
# Use a relatively stable period (inclusive) over which to take the weekday averages.
weekadjdates=[datetoday('2021-04-03'),datetoday('2021-05-14')]
weekadj=np.zeros(7)
for (date,n) in zip(apicases['date'],apicases['newCasesBySpecimenDate']):
  day=datetoday(date)
  if day>=weekadjdates[0] and day<=weekadjdates[1]: weekadj[day%7]+=n
weekadjp=weekadj*7/sum(weekadj)

if reduceltla!=None:
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
      cases[place][d]+=n/specadj[d]/weekadjp[day%7]
  places=sorted(list(cases))
else:
  places=sorted(list(vocnum))

# Restrict to places for which there is at least some of each variant, and bundle the remaining locations together as "Other"
okplaces=set([place for place in places if vocnum[place][:,0].sum()>0 and vocnum[place][:,1].sum()>0])
#okplaces=set(places)
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

if 0:
  vdir='tempd'
  for place in places:
    with open(os.path.join(vdir,place.replace(' ','_')),'w') as fp:
      for w in range(nweeks):
        print(daytodate(firstweek+voclen*w),"%6d %6d"%tuple(vocnum[place][w]),file=fp)

# Work out pre-Delta case counts, amalgamated to at least region level
if reduceltla!=None:
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
# dh = 1 standard deviation
def Rdesc(h0,dh):
  (Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-zconf*dh,h0,h0+zconf*dh]]
  return "%.0f%% (%.0f%% - %.0f%%)"%(T,Tmin,Tmax)

def Ddesc(h0,dh):
  (Dmin,D,Dmax)=[log(2)/h for h in [h0+zconf*dh,h0,h0-zconf*dh]]
  return "%.2f (%.2f - %.2f)"%(D,Dmin,Dmax)

def Gdesc(h0,dh):
  return "%.2f%% (%.2f%% - %.2f%%)"%(h0*100,(h0-zconf*dh)*100,(h0+zconf*dh)*100)

def GDdesc(h0,hlow,hhigh):
  hh=(h0,hlow,hhigh)
  s="%.2f (%.2f - %.2f) per day"%hh
  if abs(h0)>=0.01:
    if h0>0: s+=" [doubling time: %.2f (%.2f - %.2f) days]"%tuple(log(2)/h for h in [h0,hhigh,hlow])
    elif h0<0: s+=" [halving time: %.2f (%.2f - %.2f) days]"%tuple(-log(2)/h for h in hh)
  return s

print("Estimating competitive advantage using variant counts only (not case counts)")
print("============================================================================")
print()

mincount=1
l=[]
eps=1e-20
for w in range(nweeks-1):
  if len(places)<30: print(daytodate(firstweek+w*voclen),end='   ')
  for place in places:
    vn=vocnum[place]
    if (vn[w:w+2,:]>=mincount).all():
      g=((vn[w+1][1]+eps)/(vn[w+1][0]+eps))/((vn[w][1]+eps)/(vn[w][0]+eps))
      l.append(g)
      if len(places)<30:
        print(" %6.1f%%"%(log(g)/voclen*100),end='')
    else:
      if len(places)<30:
        print("  ------",end='')
  if len(places)<30: print()
l.sort()
n=len(l)
k=int(binom.ppf((1-conf)/2,n,0.5))
med=log((l[n//2]+l[(n-1)//2])/2)/voclen
low=log(l[k-1])/voclen
high=log(l[n-1-k])/voclen
print("Separate location & weeks, unweighted high-low non-parametric test: %.2f%% (%.2f%% - %.2f%%)"%(med*100,low*100,high*100))
print()

l=[]
for w in range(nweeks-1):
  for place in places:
    vn=vocnum[place]
    if (vn[w:w+2,:]>0).all():
      wt=sqrt(1/(1/vn[w:w+2,:]).sum())
      g=(vn[w+1][1]/vn[w+1][0])/(vn[w][1]/vn[w][0])
      l.append((g,wt))
l.sort()
wts=np.array([wt for (g,wt) in l])
n=len(l)
nsamp=int(1e6/len(places))
rand=bernoulli.rvs(0.5,size=[nsamp,n])
samp=rand@wts
samp.sort()
wtlow=samp[int(nsamp*(1-conf)/2)]
wtmed=samp[int(nsamp/2)]
wthigh=samp[int(nsamp*(1+conf)/2)]
wt=0
low=med=high=None
for i in range(n):
  wt+=l[i][1]
  if low==None and wt>wtlow-1e-6: low=log(l[i][0])/voclen
  if med==None and wt>wtmed-1e-6: med=log(l[i][0])/voclen
  if high==None and wt>wthigh-1e-6: high=log(l[i][0])/voclen
print("Separate location & weeks, weighted high-low non-parametric test: %.2f%% (%.2f%% - %.2f%%)"%(med*100,low*100,high*100))
print()

for w in range(nweeks-1):
  l=[]
  for place in places:
    vn=vocnum[place]
    if (vn[w:w+2,:]>0).all():
      wt=sqrt(1/(1/vn[w:w+2,:]).sum())
      g=(vn[w+1][1]/vn[w+1][0])/(vn[w][1]/vn[w][0])
      l.append((g,wt))
  if l==[]: continue
  l.sort()
  wts=np.array([wt for (g,wt) in l])
  n=len(l)
  nsamp=int(1e6/len(places))
  rand=bernoulli.rvs(0.5,size=[nsamp,n])
  samp=rand@wts
  samp.sort()
  wtlow=samp[int(nsamp*(1-conf)/2)]
  wtmed=samp[int(nsamp/2)]
  wthigh=samp[int(nsamp*(1+conf)/2)]
  wt=0
  low=med=high=None
  for i in range(n):
    wt+=l[i][1]
    if low==None and wt>wtlow-1e-6: low=log(l[i][0])/voclen
    if med==None and wt>wtmed-1e-6: med=log(l[i][0])/voclen
    if high==None and wt>wthigh-1e-6: high=log(l[i][0])/voclen
  print(daytodate(firstweek+voclen*w),"Separate locations, weighted high-low non-parametric test: %6.2f%% (%6.2f%% - %6.2f%%)"%(med*100,low*100,high*100))
print()

from scipy.special import betaln
from scipy.integrate import quad
from scipy import inf
def crossratiosubdivide(matgen,duration=voclen):
  tot=np.zeros([2,2],dtype=int)
  ndiv=20
  logp=np.zeros(ndiv)
  L0=L1=0
  for M in matgen:
    tot+=M
    if (M>0).all():
      c=1/((1/M.flatten()).sum())
      T=M[0,0]*M[1,1]/(M[0,1]*M[1,0])
      L0+=c*log(T);L1+=c
      for i in range(ndiv):
        x=(hmin+(i+.5)/ndiv*(hmax-hmin))*duration# Convert to weekly growth rate
        a,b,c,d=M[0,0],M[0,1],M[1,0],M[1,1]
        l0=d*x-(betaln(a,b)+betaln(c,d))
        # Faff around finding maximum to avoid underflow in integral
        e=exp(x)
        X=b+d-1;Y=a+b;Z=c+d
        A=Y+Z-X
        B=Y/e+Z-X*(1+1/e)
        C=-X/e
        z0=(-B+sqrt(B**2-4*A*C))/(2*A)
        l1=(b+d-1)*log(z0) - (a+b)*log(1+z0) - (c+d)*log(1+e*z0)
        sc=sqrt((b+d-1)/z0**2-(a+b)/(1+z0)**2-(c+d)*e**2/(1+e*z0)**2)# Set inverse length scale
        #res=quad(lambda z: exp( (b+d-1)*log(z) - (a+b)*log(1+z) - (c+d)*log(1+e*z) - l1 ), 0, inf)
        res=quad(lambda x: exp( (b+d-1)*log(x/sc) - (a+b)*log(1+x/sc) - (c+d)*log(1+e*x/sc) - l1 ), 0, inf)
        logp[i]+=log(res[0]/sc)+l0+l1
  if (tot==0).any():
    print("Can't estimate best transmission factor because VOC count matrix has 1 or more zero entries");return
  g=log(tot[0,0]*tot[1,1]/(tot[0,1]*tot[1,0]))/duration
  dg=sqrt((1/tot.flatten()).sum())/duration
  print("Overall cross ratio:",Gdesc(g,dg),Rdesc(g,dg),tot.flatten())
  g=L0/L1/duration
  dg=1/sqrt(L1)/duration
  print("Inverse variance weighting method using log(CR):",Gdesc(g,dg),Rdesc(g,dg))
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
  print("Likelihood method using log(CR):",Gdesc(g0,dg),Rdesc(g0,dg))
  print()

# Simple regression with 1/(1/v0+1/v1) weighting
def simpleregress(NV):
  DT=np.array(range(firstweek,firstweek+voclen*nweeks,voclen))
  W=NV[:,0]*NV[:,1]/(NV.sum(axis=1)+1e-20)
  day0=DT.sum()/nweeks
  if (W>0).sum()<=2: return (day0,0,1)
  X=DT-day0
  Y=np.log((NV[:,1]+1e-20)/(NV[:,0]+1e-20))
  m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
  r=np.array([sum(W*Y),sum(W*X*Y)])
  c=np.linalg.solve(m,r)
  mi=np.linalg.pinv(m)
  R=c[0]+c[1]*X-Y
  overdis=(R*R*W).sum()/len(R)
  dg=sqrt(mi[1,1]*overdis)
  return (day0-c[0]/c[1],c[1],dg)

n=1+len(places)
def NLL_vonly(xx):
  g=xx[0]
  LL=0
  for (i,place) in enumerate(places):
    nn=vocnum[place]
    t0=xx[1+i]
    for w in range(nweeks):
      G=(w*voclen-t0)*g
      LL+=-(nn[w][0]+nn[w][1])*log(1+exp(G))+nn[w][1]*G
  return -LL/1000

if voclen>=7:
  for w in range(nweeks-1):
    day0=lastweek-(nweeks-w)*voclen+1
    print(daytodate(day0),"-",daytodate(day0+2*voclen-1))
    crossratiosubdivide(vocnum[place][w:w+2] for place in places)

print("All week pairs:")
crossratiosubdivide(vocnum[place][w:w+2] for place in places for w in range(nweeks-1))
print("First week to last week:")
crossratiosubdivide((vocnum[place][0:nweeks:nweeks-1] for place in places), duration=voclen*(nweeks-1))

if len(places)<50:
  print("--- Inverse variance weighted regression using combined counts for %s ---"%areacovered)
  sr=simpleregress(sum(vocnum.values()))
  xx=[sr[1]]# Overall growth
  print("%s:"%areacovered,Gdesc(sr[1],sr[2]),Rdesc(sr[1],sr[2]),"   crossover on",daytodate(sr[0]))
  print()
  print("--- Inverse variance weighted regression using counts for each %s ---"%locationsize)
  s0=s1=0
  for place in places:
    sr=simpleregress(vocnum[place])
    xx.append(sr[0]-firstweek)# Intercept of [place]
    print("%25s:"%place,Gdesc(sr[1],sr[2]),Rdesc(sr[1],sr[2]),"   crossover on",daytodate(sr[0]))
    iv=1/sr[2]**2
    s0+=iv
    s1+=iv*sr[1]
  print()
  g=s1/s0;dg=sqrt(1/s0)
  print("Inverse variance weighted count for %s controlling for %s:"%(areacovered,locationsize),Gdesc(g,dg),Rdesc(g,dg))
  print()
  
  bounds=[(hmin,hmax)]+[(x-100,x+100) for x in xx[1:]]
  res=minimize(NLL_vonly,xx,bounds=bounds,method="SLSQP",options=minopts)
  print("--- Quasi-Poisson regression controlling for %s ---"%locationsize)
  print("Growth: %.2f%%    crossover on"%(res.x[0]*100),daytodate(firstweek+res.x[1]))
  print()

print()
if reduceltla==None: sys.exit(0)

print("Estimating competitive advantage using variant counts together with case counts")
print("===============================================================================")
print()

# L=ndays+2*bmsig
# i is time from start, in days
# t=i/L
# growth[i] = bmscale*sqrt(L)*(t*X_0 + sum_{n=1}^{N-1} sqrt(2)/pi*exp(-(n*bmsig/L)^2/2)*sin(n*pi*t)/n*X_n)
# where X_n ~ N(0,1),  n=0,...,N; N=ceil(4*L/bmsig), say

bmL=ndays+int(bmsig*2+0.999)# Add on bmsig*2 to eliminate periodicity effects
bmN=int(2.5*bmL/bmsig+1)
bmsin=[sin(r*pi/bmL) for r in range(2*bmL)]
bmweight=[0]+[sqrt(2)/pi*exp(-(n*bmsig/bmL)**2/2)/n for n in range(1,bmN)]
bmsin2=[np.array([bmweight[n]*bmsin[(i*n)%(2*bmL)] for n in range(bmN)]) for i in range(ndays)]

# Need to scale the variables being optimised over to keep SLSQP happy
#condition=np.array([50,50,1000,1000]+[1.]*bmN)
t=40000*np.array(bmweight[1:])*bmscale/np.arange(1,bmN)+1
condition=np.array([70,80,5000,5000,t[0]]+list(t))

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
  w=bmscale*sqrt(bmL)
  GG=[]
  H=exp(h)
  for i in range(ndays-1):
    t=i/bmL
    gu=t*xx[4]+np.dot(bmsin2[i],xx[4:])
    #for n in range(1,bmN):
    #  gu+=bmweight[n]*bmsin[(i*n)%(2*bmL)]*xx[4+n]
    g=g0+w*gu
    GG.append(g)
    G=exp(g)
    AA.append(AA[-1]*G)
    BB.append(BB[-1]*G*H)
  return AA,BB,GG


# Return negative log likelihood (negative because scipy can only minimise, not maximise)
# If const is true then add in all the constant terms (that don't affect the optimisation)
def NLL(xx_conditioned,lcases,lvocnum,sig0,asc,lprecases,const=False):
  xx=xx_conditioned/condition
  tot=0
  
  # Prior on starting number of cases of non-Delta: assume starts off similar to total number of cases
  a0=log(lcases[0]+.5)
  tot+=-((xx[0]-a0)*isd0)**2/2
  if const: tot-=log(2*pi/isd0**2)/2
  
  # Very weak prior on starting number of cases of Delta
  tot+=-((xx[1]-(a0-4))*isd1)**2/2
  if const: tot-=log(2*pi/isd1**2)/2
  
  # Prior on h
  tot+=-(xx[2]*isd2)**2/2
  if const: tot-=log(2*pi/isd2**2)/2
  
  a,b=lprecases[0]+.5,lprecases[1]+.5
  g0=log(b/a)/7
  v0=(1/a+1/b)/49+sig0**2
  tot+=-(xx[3]-g0)**2/(2*v0)
  if const: tot-=log(2*pi*v0)/2
  
  AA,BB,GG=expand(xx)
  # Component of likelihood due to number of confirmed cases seen
  for i in range(ndays):
    mu=asc*(AA[i]+BB[i])
    r=mu*nif1a[i]/(1-nif1a[i])
    n=lcases[i]
    # n ~ Negative binomial(mean=mu, variance=mu/nif1)
    # max with -10000 because the expression is unbounded below which can cause a problem for SLSQP
    if model=="scaledpoisson":
      tot+=max((-mu+n*log(nif1a[i]*mu))*nif1a[i],-10000)
      if const: tot+=log(nif1a[i])-gammaln(nif1a[i]*n+1)# Approx normalisation
    elif model=="NBBB":
      tot+=max(gammaln(n+r)+r*lognif1a[i]+n*log1mnif1a[i]-gammaln(r),-10000)
      if const: tot+=-gammaln(n+1)
    elif model=="NBBB+magicprior":
      tot+=max(gammaln(n+r)-nif1a[i]*gammaln(mu+r)+n*log1mnif1a[i],-10000)
      if const: tot+=-gammaln(n+1)
    else: raise RuntimeError("Unrecognised model "+model)
  
  # Term to regulate change in growth rate
  for i in range(bmN):
    tot+=-xx[4+i]**2/2
  if const: tot-=bmN*log(2*pi)/2
  
  # Term to align the variant numbers with VOC count data
  for w in range(nweeks):
    endweek=lastweek-(nweeks-1-w)*voclen-minday
    A=sum(AA[endweek-(voclen-1):endweek+1])
    B=sum(BB[endweek-(voclen-1):endweek+1])
    f=nif2/(1-nif2);a=f*A;b=f*B
    r,s=lvocnum[w][0],lvocnum[w][1]
    if model=="scaledpoisson":
      r1,s1=nif2*r,nif2*s
      tot+=r1*log(A/(A+B))+s1*log(B/(A+B))
      if const: tot+=gammaln(r1+s1+1)-gammaln(r1+1)-gammaln(s1+1)+log(nif2)
    elif model=="NBBB" or model=="NBBB+magicprior":
      if abs(a+b)<10000*(r+s):
        tot+=gammaln(a+r)+gammaln(b+s)-gammaln(a+b+r+s)+gammaln(a+b)-gammaln(a)-gammaln(b)
      else:
        tot+=r*log(A/(A+B))+s*log(B/(A+B))
      if const: tot+=gammaln(r+s+1)-gammaln(r+1)-gammaln(s+1)
    else: raise RuntimeError("Unrecognised model "+model)

  return -tot

def Hessian(xx,lcases,lvocnum,sig0,asc,lprecases):
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
        v+=s1*s2*NLL(x*condition,lcases,lvocnum,sig0,asc,lprecases)
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
  return H
      
# Returns log likelihood
def optimiseplace(place,hint=np.zeros(bmN+4),fixedh=None,statphase=False):
  xx=np.copy(hint)
  # bounds[2][0]=0 prejudges Delta as being at least as transmissible as B.1.1.7. This helps SLSQP not get stuck in some cases
  # though would need to relax this constraint if dealing with other variants where it might not be true.
  bounds=[(-10,20),(-10,20),(0,1),(-1,1)]+[(-10,10)]*bmN
  if fixedh!=None: xx[2]=fixedh;bounds[2]=(fixedh,fixedh)
  res=minimize(NLL,xx*condition,args=(cases[place],vocnum[place],sig0,asc,precases[prereduce(place)]),bounds=bounds*np.repeat(condition,2).reshape([len(bounds),2]),method="SLSQP",options=minopts)
  if not res.success:
    print(res)
    print(fixedh)
    print(place)
    print("xx =",xx)
    print("bounds =",bounds)
    print("lcases =",list(cases[place]))
    print("lprecases =",precases[prereduce(place)])
    print("lvocnum =",vocnum[place])
    for x in sorted(list(opts)): print("%s:"%x,opts[x])
    print("nweeks, ndays, minday, lastweek =",nweeks,",",ndays,",",minday,",",lastweek)
    raise RuntimeError(res.message)
  xx=res.x/condition

  # Work out log likelihood including constant terms
  LL=-NLL(res.x,cases[place],vocnum[place],sig0,asc,precases[prereduce(place)],const=True)
  
  # If 'statphase', make the log likelihood a better approximation to log(integral over all parameters) using stationary phase approximation
  if statphase:
    H=Hessian(xx,cases[place],vocnum[place],sig0,asc,precases[prereduce(place)])
    det=np.linalg.det(H)
    N=H.shape[0]
    if det<=0: print("Warning: Hessian not positive for %s. Can't make corrected log likelihood."%place);det=1
    LL+=N*log(2*pi)/2-log(det)/2

  # Return optimum xx log likelihood
  return res.x/condition,LL

def getsamples(place,xx0):
  H=Hessian(xx0,cases[place],vocnum[place],sig0,asc,precases[prereduce(place)])
  Hcond=H/condition/condition[:,None]
  eig=np.linalg.eigh(Hcond)
  # np.diag(np.matmul(np.matmul(np.transpose(eig[1]),Hcond),eig[1])) ~= eig[0]
  if not (eig[0]>0).all(): print("Hessian not +ve definite so can't do full confidence calculation");return None
  nsamp=10000
  N=bmN+4
  t=norm.rvs(size=[nsamp,N])# nsamp x N
  sd=eig[0]**(-.5)# N
  u=t*sd# nsamp x N
  samp_cond=np.matmul(u,np.transpose(eig[1]))# nsamp x N
  samp=samp_cond/condition
  SSS=np.zeros([nsamp,2,ndays])
  for i in range(nsamp):
    xx=xx0+samp[i]
    AA,BB,GG=expand(xx)
    SSS[i,0,:]=AA
    SSS[i,1,:]=BB
  return SSS


# Generate sample paths conditional on h = xx0[2] + dhsamp
# xx0 = N-vector that is optimal conditioned on xx0[2]
# dhsamp = nsamp-vector of deltas in xx0[2] to sample over (an externally imposed normal - not from the log likelihood)
# Use normal approximation with Hessian from log likelihood for other co-ordinates
#
# Types in terms of dimensions and conditionedness:
#                       L^this  condition^this
# xx, xx0, dhsamp, samp      1      0
# condition                  0      1
# H                         -2      0
# Hcond, Hcond_, Hcond__    -2     -2
# eig[0]                    -2     -2
# eig[1], t, s1              0      0
# sd, u, s0                  1      1
def getcondsamples(place,xx0,dhsamp):
  N=bmN+4
  H=Hessian(xx0,cases[place],vocnum[place],sig0,asc,precases[prereduce(place)])
  Hcond=H/condition/condition[:,None]
  Hcond__=np.delete(np.delete(Hcond,2,0),2,1)# N-1 x N-1
  Hcond_=np.delete(Hcond[2],2,0)# N-1
  eig=np.linalg.eigh(Hcond__)
  # np.matmul(np.matmul(np.transpose(eig[1]),Hcond__),eig[1]) ~= np.diag(eig[0])
  m=eig[0].min()
  if m<=0: print("Hessian not +ve definite in getcondsamples so can't do full confidence calculation");return None
  if m<1e-6: print("Warning: Hessian has a very low eigenvalue in getcondsamples:",m)
  nsamp=len(dhsamp)
  t=norm.rvs(size=[nsamp,N-1])# nsamp x N-1
  sd=eig[0]**(-.5)# N-1
  u=t*sd# nsamp x N-1
  s0=np.insert(np.matmul(u,np.transpose(eig[1])),2,0,1)# nsamp x N
  s1=np.insert(-np.linalg.solve(Hcond__,Hcond_),2,1,0)# N
  samp=s0/condition+dhsamp[:,None]*s1# nsamp x N
  SSS=np.zeros([nsamp,2,ndays])
  # This is the slow bit.
  for i in range(nsamp):
    xx=xx0+samp[i]
    AA,BB,GG=expand(xx)
    SSS[i,0,:]=AA
    SSS[i,1,:]=BB
  return SSS
  
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

def fullprint(AA,BB,lvocnum,lcases,T=None,Tmin=None,Tmax=None,area=None,using='',samples=None):
  print("ModV1    = modelled number of new cases of Alpha on this day multiplied by the ascertainment rate")
  print("ModV2    = modelled number of new cases of Delta on this day multiplied by the ascertainment rate")
  print("Pred     = predicted number of cases seen this day = ModV1+ModV2")
  print("Seen     = number of cases observed this day, after weekday adjustment, from api/dashboard")
  print("PredV1   = p*Pred, where p = proportion of Alpha amongst obserbed variant counts from",source)
  print("PredV2   = (1-p)*Pred")
  print("SeenV1   = p*Seen")
  print("SeenV2   = (1-p)*Seen")
  print("Q        = estimated reproduction rate of Alpha on this day")
  print("R        = estimated reproduction rate of Delta on this day")
  print("ModV1min = ModV1 min confidence interval")
  print("ModV1med = ModV1 mode confidence interval")
  print("ModV1max = ModV1 max confidence interval")
  print("ModV2min = ModV1 min confidence interval")
  print("ModV2med = ModV1 mode confidence interval")
  print("ModV2max = ModV1 max confidence interval")
  print("Qmin     = Q min confidence interval")
  print("Qmed     = Q mode confidence interval")
  print("Qmax     = Q max confidence interval")
  print("Rmin     = R min confidence interval")
  print("Rmed     = R mode confidence interval")
  print("Rmax     = R max confidence interval")
  print()
  ave=1# Number of days over which to average to get growth rate estimates
  # samples[ sample number, 0 or 1 for alpha/delta, day ] = cases
  nsamp=samples.shape[0]
  nmin=int((1-conf)/2*nsamp)
  nmed=nsamp//2
  nmax=int((1+conf)/2*nsamp)
  sa=np.sort(samples,axis=0)# For each variant, day, put samples in order of cases
  sr0=np.log(samples[:,:,ave:]/samples[:,:,:-ave])/ave# sr, sr0 = growths g(V1), g(V2)
  sr=np.sort(sr0,axis=0)# For each variant, day, put samples in order of scaled change from one day to next (R-number)
  tr=np.sort(sr0[:,1,:]-sr0[:,0,:],axis=0)# For each day, d, put samples in order of g(newvar,d)-g(oldvar,d)
  QQ=sr[:,0,-1]
  RR=sr[:,1,-1]
  TT=list(tr[:,-1])
  if area!=None and args.graph_filename!=None:
    graphdata=sanitise(args.graph_filename+'_'+area+'.dat')
    graphfp=open(graphdata,'w')
  else:
    graphfp=None
  def mprint(*a,**b):
    print(*a,**b)
    if graphfp!=None: print(*a,**b,file=graphfp)
  mprint("#     Date     ModV1     ModV2      Pred      Seen      PredV1    PredV2    SeenV1    SeenV2          Q       R  ModV1min  ModV1med  ModV1max  ModV2min  ModV2med  ModV2max    Qmin    Qmed    Qmax    Rmin    Rmed    Rmax")
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
    #mprint(" %12g %12g"%(asc*(sa[nmed][0][i]-AA[i]),asc*(sa[nmed][1][i]-BB[i])),end='')
    if i<ndays-ave:
      Q,R=(log(AA[i+ave]/AA[i])/ave,log(BB[i+ave]/BB[i])/ave)
      mprint("   %7.4f %7.4f"%(Q,R),end='')
    else:
      mprint("         -       -",end='')
    mprint(" %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f"%(asc*sa[nmin,0,i],asc*sa[nmed,0,i],asc*sa[nmax,0,i],asc*sa[nmin,1,i],asc*sa[nmed,1,i],asc*sa[nmax,1,i]),end='')
    if i<ndays-ave:
      mprint(" %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f"%(sr[nmin,0,i],sr[nmed,0,i],sr[nmax,0,i],sr[nmin,1,i],sr[nmed,1,i],sr[nmax,1,i]))
    else:
      mprint("       -       -       -       -       -       -")
  EQ="Estd growth in %s: "%nonvariant+GDdesc(QQ[nmed],QQ[nmin],QQ[nmax])
  ER="Estd growth in %s: "%variant+GDdesc(RR[nmed],RR[nmin],RR[nmax])
  # Note that T is not 100(R/Q-1) here because AA, BB are derived from a sum of locations each of which has extra transm T,
  # but because of Simpson's paradox, that doesn't mean the cross ratio of AAs and BBs is also T.
  ETA="Estd growth advantage = "
  if T!=None: ETA+=GDdesc(T,Tmin,Tmax)
  else: ETA+=GDdesc(TT[nmed],TT[nmin],TT[nmax])
  if area!=None: print("Summary");print(area+using)
  print(EQ)
  print(ER)
  print(ETA)
  ETA+="\\n(CIs show within-model statistical uncertainty; model assumptions lead to other uncertainties)"
  print()
  if graphfp!=None:
    graphfp.close()
    now=datetime.utcnow().strftime('%Y-%m-%d')
    po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
    # Use this write function to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
    def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))
    
    write('set terminal pngcairo font "sans,13" size 1920,1280')
    write('set bmargin 7;set lmargin 13;set rmargin 13;set tmargin 5')
    write('set xdata time')
    write('set timefmt "%Y-%m-%d"')
    write('set format x "%Y-%m-%d"')
    write('set xtics nomirror rotate by 45 right offset 0.5,0')
    write('set label "Location: %s\\nAs of %s:\\n%s\\n%s\\n%s" at screen 0.45,0.9'%(area+using,daytodate(minday+ndays-ave-1),EQ,ER,ETA))

    mi=[1e9]*50;mx=[-1e9]*50
    with open(graphdata,'r') as fp:
      for x in fp:
        if x[0]!='#':
          y=x.split()
          for (i,z) in enumerate(y):
            try: f=float(z);mi[i]=min(mi[i],f);mx[i]=max(mx[i],f)
            except: pass
    
    graphfnR=sanitise(args.graph_filename+'_'+area+'_R.png')
    mi0=min(mi[17],mi[18],mi[20],mi[21])
    mx0=max(mx[17],mx[18],mx[20],mx[21])
    yscale='[:] [:%f]'%(mi0+(mx0-mi0)/.8)
    write('set output "%s"'%graphfnR)
    write('set key left')
    write('set ylabel "Estimated growth rate')
    write('set style fill transparent solid 0.25')
    write('set style fill noborder')
    #write('set y2tics')
    write('set title "Estimated continuous growth rates of %s and %s variants in %s\\n'%(nonvariant,variant,area+using)+
          'Fit made on %s using https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/vocfit.py\\n'%now+
          'Data sources: %s; https://coronavirus.data.gov.uk/"'%fullsource)
    write('plot %s "%s" u 1:19 w lines lc "green" lw 3 title "Estimated growth in %s", "%s" u 1:18:20 with filledcurves lc "green" title "", "%s" u 1:22 w lines lc "blue" lw 3 title "Estimated growth in %s", "%s" u 1:21:23 with filledcurves lc "blue" title ""'%(yscale,graphdata,nonvariant,graphdata,graphdata,variant,graphdata))
    print("Written graph to %s"%graphfnR)
    
    for yaxis in ["lin","log"]:
      graphfn=sanitise(args.graph_filename+'_'+area+'_'+yaxis+'.png')
      if yaxis=="log": write('set logscale y');write('unset y2tics');yscale='[:] [0.99:%f]'%(max(mx[1:7])**(1/.8))
      else: yscale='[:] [0:%f]'%(max(mx[1:7])/.8)
      write('set key top left')
      write('set output "%s"'%graphfn)
      write('set ylabel "New cases per day (scaled down to match ascertainment rate of %0.f%%)"'%(100*asc))
      write('set title "Estimated new cases per day of %s and %s in %s\\n'%(nonvariant,variant,area+using)+
            'Fit made on %s using https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/vocfit.py\\n'%now+
            'Data sources: %s; https://coronavirus.data.gov.uk/"'%fullsource)
      write('plot %s "%s" u 1:2 with lines lw 3 title "Modelled %s", "%s" u 1:3 with lines lw 3 title "Modelled %s", "%s" u 1:4 with lines lw 3 title "Modelled total", "%s" u 1:5 with lines lt 6 lw 3 title "Confirmed cases (all variants, weekday adjustment)", "%s" u 1:6 lt 1 pt 6 lw 3 title "Proportion of %s scaled up to modelled total", "%s" u 1:7 lt 2 pt 6 lw 3 title "Proportion of %s scaled up to modelled total"'%(yscale,graphdata,nonvariant,graphdata,variant,graphdata,graphdata,graphdata,nonvariant,graphdata,variant))
      print("Written graph to %s"%graphfn)
    
    p.close();po.wait()
  return Q,R

def printsummary(summary):
  print("Location                       Q     R      T")
  for place in places:
    (Q,R,T,Tmin,Tmax)=summary[place]
    print("%-25s  %5.2f %5.2f  %4.0f%%"%(place,Q,R,T),end='')
    if Tmin!=None: print(" ( %4.0f%% - %4.0f%% )"%(Tmin,Tmax))
    else: print()
  print()
  print("Q = point estimate of reproduction rate of %s on"%nonvariant,daytodate(maxday-1))
  print("R = point estimate of reproduction rate of %s on"%variant,daytodate(maxday-1))
  print("T = estimated competitive advantage = R/Q as a percentage increase")
  print()

if mode=="local growth rates":
  summary={}
  for place in places:
    printplaceinfo(place)
    xx0,L0=optimiseplace(place)
    if len(places)<10:
      xx1,L1=optimiseplace(place,statphase=True)
      print("Corrected log likelihood",L1)
    AA,BB,GG=expand(xx0)
    h0=xx0[2]
    ff=[0,L0,0]
    eps=0.01
    for i in [-1,1]:
      xx,L=optimiseplace(place,hint=xx0,fixedh=h0+i*eps)
      ff[i+1]=L
    # Use observed Fisher information to make confidence interval
    fi=(-ff[0]+2*ff[1]-ff[2])/eps**2
    if fi>0:
      dh=1/sqrt(fi)
      (Tmin,T,Tmax)=[h0-zconf*dh,h0,h0+zconf*dh]
    else:
      (Tmin,T,Tmax)=[None,h0,None]
    SSS=getsamples(place,xx0)
    print("Locally optimised growth advantage")
    Q,R=fullprint(AA,BB,vocnum[place],cases[place],area=ltla2name.get(place,place),samples=SSS)
    summary[place]=(Q,R,T,Tmin,Tmax)
  print()
  printsummary(summary)

if (type(mode)==tuple or type(mode)==list) and mode[0]=="fixed growth rate":
  summary={}
  for place in places:
    printplaceinfo(place)
    h0=mode[1]
    xx0,L0=optimiseplace(place,fixedh=h0)
    AA,BB,GG=expand(xx0)
    T=h0
    print("Predetermined growth advantage")
    SSS=getcondsamples(place,xx0,[0]*10000)
    Q,R=fullprint(AA,BB,vocnum[place],cases[place],area=ltla2name.get(place,place),samples=SSS)
    summary[place]=(Q,R,T,None,None)
  print()
  printsummary(summary)
  
if mode=="global growth rate":
  ndiv=11
  logp=np.zeros(ndiv)
  for place in places:
    printplaceinfo(place)
    print("    h     T    log lik")
    xx,L0=optimiseplace(place)
    for i in range(ndiv):
      h=(hmin+(hmax-hmin)*i/(ndiv-1))
      xx,L=optimiseplace(place,hint=xx,fixedh=h)
      logp[i]+=L-L0
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
  irange=1/sqrt(c)
  h0=(hmin+(hmax-hmin)*imax/(ndiv-1))
  dh=(hmax-hmin)*irange/(ndiv-1)
  #(Tmin,T,Tmax)=[(exp(h*mgt)-1)*100 for h in [h0-zconf*dh,h0,h0+zconf*dh]]
  (Tmin,T,Tmax)=[h0-zconf*dh,h0,h0+zconf*dh]
  print("Combined growth advantage per day: "+GDdesc(T,Tmin,Tmax))
  print()
  
  print("Re-running using global optimum growth advantage",h0)
  print()
  summary={}

  # Combine places into a single top level 'areacovered', and also possibly into regions
  # Make a set of [name of aggregate location, set of ltlas that are in it]
  regions=set(ltla2region.values())
  makeregions=(locationsize=="LTLA" and areacovered not in regions)*0# temporarily disable
  if makeregions: combinedplaces=sorted(list(regions))
  else: combinedplaces=[]
  combinedplaces.append(areacovered)
  TSS0={};TSSS={}
  for loc in combinedplaces:
    TSS0[loc]=np.zeros([2,ndays])
    TSSS[loc]=np.zeros([nsamp,2,ndays])
  
  dhsamp=dh*norm.rvs(size=nsamp)
  TLL=0
  for place in places:
    using=' (using information from '+ltla2name.get(areacovered,areacovered)+')'
    printplaceinfo(place,using=using)
    xx0,L0=optimiseplace(place,fixedh=h0,statphase=True)
    TLL+=L0
    AA0,BB0,GG0=expand(xx0)
    TSS0[areacovered][0,:]+=AA0;TSS0[areacovered][1,:]+=BB0
    if makeregions: reg=ltla2region[place];TSS0[0,reg]+=AA0;TSS0[1,reg]+=BB0
    SSS=getcondsamples(place,xx0,dhsamp)
    if not SSS is None:
      assert SSS[:,1,-2:].max()<1e20
      TSSS[areacovered]+=SSS
      if makeregions: TSSS[reg]+=SSS
    print("Globally optimised growth advantage")
    area=None
    if place!=areacovered and (locationsize!="LTLA" or place in specialinterest or len(places)<50): area=ltla2name.get(place,place)
    Q,R=fullprint(AA0,BB0,vocnum[place],cases[place],area=area,using=using,samples=SSS)
    summary[place]=(Q,R,T,None,None)
  print()
  printsummary(summary)
  print("Corrected log likelihood",TLL)
  
  print("Total predicted counts using global optimum growth advantage")
  print()

  for loc in combinedplaces:
    print("Combined results for %s using globally optimised growth advantage"%loc)
    Q,R=fullprint(TSS0[loc][0,:],TSS0[loc][1,:],sum(vocnum.values()),[sum(cases[place][i] for place in places) for i in range(ndays)],T,Tmin,Tmax,area=ltla2name.get(loc,loc),samples=TSSS[loc])
  
  print("Combined growth advantage per day: "+GDdesc(h0,h0-zconf*dh,h0+zconf*dh))
