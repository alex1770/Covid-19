from stuff import *
import sys,re,argparse,pickle,csv
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.special import gammaln
from math import log,exp,sqrt,sin,pi
import numpy as np
from subprocess import Popen,PIPE
from datetime import datetime
import pytz

# ltla_age.csv from https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesBySpecimenDateAgeDemographics&format=csv
# wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=vaccinationsAgeDemographics&format=csv' -O vaccinedata.csv
# Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv

# Case age bands (all English LTLAs): 0-5, 5-10, ..., 85-90, 90+ (18 5-year bands and 90+)
# Vax and NIMS-population age bands for English LTLAs: 18-25, 25-30, 30-35, ..., 85-90, 90+ (18-25, 13 5-year bands, 90+)
# Standardise on these age bands: 0-5, 5-10, 10-15, 15-18, 18-25, 25-30, ..., 85-90, 90+
# by (i) zero-extending NIMS data, and
#    (ii) using simple interpolation to convert case data from 15-20, 20-25 to 15-18, 18-25
# These bands are hard-coded, perforce, including the assumption that there are equal number of case age bands and unified age bands (17)

# Unified age bands:
num2age=['00_04','05_09','10_14','15_17','18_24']+['%02d_%02d'%(i,i+4) for i in range(25,80,5)]+['80+']
numage=len(num2age)
# Basic agestring -> ageband
age2num={a:i for (i,a) in enumerate(num2age)}
# Add decodes for vax ages
for (a,age) in enumerate(num2age):
  if '_' in age: age2num[age.replace('_','-')]=a
# Flag unwanted age bands
age2num['00_59']=None
age2num['60+']=None
age2num['unassigned']=None
# Fuse case 80-84, 85-89, 90+ into 80+
age2num['80_84']=numage-1
age2num['85_89']=numage-1
age2num['90+']=numage-1
# Temporarily assign NIMS/ONS population U18s to band 0
age2num['Under 18']=0
# Temporarily assign case 15_19 and 20_24 to 15_17 and 18_24 respectively
age2num['15_19']=age2num['15_17']
age2num['20_24']=age2num['18_24']

# Aim for numpy arrays:
# ltlapop[ltla][age]
# vax1[day][ltla][age]  (total 1st dose vaccines by day `day`)
# vax2[day][ltla][age]  (total 2nd dose vaccines by day `day`)
# cumcases[day][ltla][age]  (total cases by day `day`)
# newcases[day][ltla][age]  (new cases on day `day`)
# sanger[] TBD

def sanitise(fn): return fn.replace(' ','_').replace("'","")

ltlaengdata=loadcsv("Local_Authority_District_to_Region_(April_2021)_Lookup_in_England.csv")
ltla2region=dict(zip(ltlaengdata['LAD21CD'],ltlaengdata['RGN21NM']))
ltla2name=dict(zip(ltlaengdata['LAD21CD'],map(sanitise,ltlaengdata['LAD21NM'])))

# 1. Sanger uses LAD19 except E06000053 (Isles of Scilly) isn't present. I assume it's fused into E06000052 (Cornwall).
# 2. Dashboard/api (including NIMS/population/vax data) uses LAD19 except that E09000001 (City of London) has been fused into E09000012 (Hackney).
#    and E06000053 (Isles of Scilly) is fused into E06000052 (Cornwall). Call this LAD19+.
# 3. NIMS from NHS page, prior to July 2021 uses LAD20 (cf LAD19: E0700000[4567] fused into E06000060).
# 4. NIMS from NHS page, from July 2021 uses LAD21 (cf LAD20: E0700015[0236] fused into E06000061, E0700015[145] fused into E06000062).
# Could almost work solely from dashboard (incl populations), in which case could standardise on LAD19+, but then wouldn't have under 18 population, and perhaps more importantly
# we wouldn't have the option to use ONS populations. So let's standardise on LAD21+ (=LAD21 with E09000001 fused into E09000012 and E06000053 fused into E06000052),
# and get the populations (either NIMS or ONS) from the NHS page. Those only go up to age 80+, so let's also fuse age groups 80-84, 85-90, 90+ into 80+.
fuseltla=dict(zip(ltlaengdata['LAD21CD'],ltlaengdata['LAD21CD']))
# api/dashboard quirks (call this LAD19+)
fuseltla['E09000001']='E09000012'
fuseltla['E06000053']='E06000052'
# fusings to use LAD20+
fuseltla['E07000004']='E06000060'
fuseltla['E07000005']='E06000060'
fuseltla['E07000006']='E06000060'
fuseltla['E07000007']='E06000060'
# fusings to use LAD21+
fuseltla['E07000150']='E06000061'
fuseltla['E07000152']='E06000061'
fuseltla['E07000153']='E06000061'
fuseltla['E07000156']='E06000061'
fuseltla['E07000151']='E06000062'
fuseltla['E07000154']='E06000062'
fuseltla['E07000155']='E06000062'
# Convert standardised (fused) ltlas to and from numbers:
num2ltla=sorted(list(set(fuseltla.values())))
ltla2num={ltla:i for (i,ltla) in enumerate(num2ltla)}
numltla=len(num2ltla)
s=set(ltlaengdata['LAD21CD']);assert all(x in s for x in fuseltla.values())

indir='inputs'
os.makedirs(indir,exist_ok=True)
apiday0=apiday()

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



### Options ###

source="Sanger"

r0=0.3

vaxeffecttime=20# Days before vaccine is presumed to have a decent effect

ltlaexclude=set()
#ltlaexclude=set(['E08000001','E12000002'])# Bolton, Manchester
ltlaset="All"
#ltlaset="London"
#ltlaset="Bolton"
#ltlaset="Hartlepool"

mgt=5# Mean generation time in days

# Earliest day to use VOC count data, given as end-of-week. Will be rounded up to match same day of week as lastweek.
firstweek=datetoday('2021-05-01')

# Case ascertainment rate
car=0.4

# Discard this many cases at the end of the list of cases by specimen day
discardcasedays=2# Make need to increase to 3 to allow for Wales and NI late reporting at weekends and bank holidays
missingcasedays=4# Expect 4 missing days when using newCasesBySpecimenDateAgeDemographics; setting missingcasedays at 4 means it will refresh cases every day
missingvaxdays=12# Only expect 0 missing days using vaccinationsAgeDemographics, but we only need to refresh every 12 days or so because vaccines don't have an immediate effect

minopts={"maxiter":10000,"eps":1e-4}

voclen=(1 if source=="COG-UK" else 7)

conf=0.95

### End options ###

opts={
  "VOC source": source,
  "LTLA set": ltlaset,
  "LTLA exclude": list(ltlaexclude),
  "Generation time (days)": mgt,
  "Earliest week (using end of week date) to use VOC count data": daytodate(firstweek),
  "Case ascertainment rate": car,
  "Number of days of case data to discard": discardcasedays,
  "Minimiser options": minopts,
  "Length of time period over which VOC counts are given (days)": voclen,
  "Confidence level": conf,
}

if args.save_options!=None:
  with open(args.save_options,'w') as fp: json.dump(opts,fp,indent=2)

if args.load_options!=None:
  with open(args.load_options,'r') as fp: lopts=json.load(fp)
  for x in lopts: opts[x]=lopts[x]

source=opts["VOC source"]
ltlaset=opts["LTLA set"]
ltlaexclude=set(opts["LTLA exclude"])
mgt=opts["Generation time (days)"]
firstweek=datetoday(opts["Earliest week (using end of week date) to use VOC count data"])
car=opts["Case ascertainment rate"]
discardcasedays=opts["Number of days of case data to discard"]
minopts=opts["Minimiser options"]
voclen=opts["Length of time period over which VOC counts are given (days)"]
conf=opts["Confidence level"]

zconf=norm.ppf((1+conf)/2)

np.set_printoptions(precision=3,linewidth=200)

print("Options:")
print()
for x in sorted(list(opts)): print("%s:"%x,opts[x])
print()
sys.stdout.flush()

okplaces=set(ltla for ltla in num2ltla if includeltla(ltla,ltlaset) and not ltla in ltlaexclude)
places=sorted(list(okplaces))
npl=len(places)
# ltla2num={ltla:i for (i,ltla) in enumerate(places)}

regions=sorted(list(set(ltla2region.values())))
nreg=len(regions)
reg2num={reg:i for (i,reg) in enumerate(regions)}
ltla2regnum=np.array([reg2num[ltla2region[ltla]] for ltla in num2ltla])

# Get VOC data
# Alter: get this to auto-update
if source=="Sanger":
  fullsource="Wellcome Sanger Institute"
  assert voclen==7
  sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
  
  lastweek=datetoday(max(sanger['WeekEndDate']))
  nweeks=(lastweek-firstweek)//voclen+1
  # Sanger week number is nweeks-1-(lastweek-day)//voclen

  # Get Sanger (variant) data into a numpy array
  vocnum=np.zeros([nweeks,numltla,2],dtype=int)
  rvocnum=np.zeros([nweeks,nreg,2],dtype=int)
  for (date,lad19,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
    day=datetoday(date)
    week=nweeks-1-(lastweek-day)//voclen
    if week>=0 and week<nweeks:
      ltla=fuseltla[lad19]
      if ltla not in okplaces: continue
      i=ltla2num[ltla]
      vocnum[week,i,int(var=="B.1.617.2")]+=n
      r=ltla2regnum[i]
      rvocnum[week,r,int(var=="B.1.617.2")]+=n
else: raise RuntimeError("Unrecognised source: "+source)

tvocnum=rvocnum.sum(axis=1)

# Get population data
ltlapopdata=loadcsv("LTLA-NIMS-populations.LAD21.18+.csv")
ltlapop=np.zeros([numltla,numage])
regionpop=np.zeros([nreg,numage])
n=len(ltlapopdata['LTLA Code'])
for age in ltlapopdata:
  a=age2num.get(age)
  if a!=None:
    for (lad21,nump) in zip(ltlapopdata['LTLA Code'],ltlapopdata[age]):
      i=ltla2num[fuseltla[lad21]]
      ltlapop[i,a]+=nump
# Estimate bands 0-5, 5-10, 10-15, 15-18 by assuming each year has an equal number
t=np.copy(ltlapop[:,0])
ltlapop[:,0]=5*t/18
ltlapop[:,1]=5*t/18
ltlapop[:,2]=5*t/18
ltlapop[:,3]=3*t/18
for i in range(numltla):
  regionpop[ltla2regnum[i],:]+=ltlapop[i,:]
rpopratio=np.array([ltlapop[i]/regionpop[ltla2regnum[i]] for i in range(numltla)])
tpopratio=ltlapop/ltlapop.sum(axis=0)

# Get case data into cases[daynum,ltlanum,agebandnum]
# The saved (pickled) format of case ages uses standardised (fused) ltlas, and standardised age bands
fn=os.path.join(indir,'ltla_age_cases.pickle')
ok=0
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    cases,caseday0=pickle.load(fp)
  ndays,nltlas,nages=cases.shape
  if caseday0+ndays>=apiday0-missingcasedays: ok=1
if not ok:
  print("Loading case-age data from api")
  la=loadcsv_it(api_v2('areaType=ltla&metric=newCasesBySpecimenDateAgeDemographics&format=csv').text.split('\n'))
  print("Processing case-age data")
  caseday0=datetoday(min(la['date']))
  caseday1=datetoday(max(la['date']))
  n=caseday1-caseday0+1
  cases=np.zeros([n,numltla,numage])
  for lad19,date,age,ncases in zip(la['areaCode'],la['date'],la['age'],la['cases']):
    a=age2num[age]
    if a!=None:
      ltla=fuseltla[lad19]
      cases[datetoday(date)-caseday0,ltla2num[ltla],a]+=ncases
  t=cases[:,:,3]+cases[:,:,4]# (15-20) + (20-25)
  cases[:,:,3]=3/10*t# Estimate of 15-18
  cases[:,:,4]=7/10*t# Estimate of 18-25
  with open(fn,'wb') as fp:
    pickle.dump((cases,caseday0),fp)

# Get vax data into vax[,,,]
# The saved (pickled) format of case ages uses standardised (fused) ltlas, and standardised age bands
fn=os.path.join(indir,'vax_age_cases.pickle')
ok=0
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    vax,vaxday0=pickle.load(fp)
  ndoses,ndays,nvaxltlas,nvaxages=vax.shape
  assert nvaxages==nages and nvaxltlas==nltlas
  if vaxday0+ndays>=apiday0-missingvaxdays: ok=1
if not ok:
  print("Loading vax-age data from api")
  la=loadcsv_it(api_v2('areaType=ltla&metric=vaccinationsAgeDemographics&format=csv').text.split('\n'))
  #la=loadcsv("vaccinedata.csv")
  print("Processing vax-age data")
  vaxday0=datetoday(min(la['date']))
  vaxday1=datetoday(max(la['date']))
  n=vaxday1-vaxday0+1
  vax=np.zeros([2,n,numltla,numage],dtype=int)
  for lad21,date,age,dose1,dose2 in zip(la['areaCode'],la['date'],la['age'],la['cumPeopleVaccinatedFirstDoseByVaccinationDate'],la['cumPeopleVaccinatedSecondDoseByVaccinationDate']):
    if lad21 in fuseltla:
      a=age2num[age]
      if a!=None:
        d=datetoday(date)-vaxday0
        i=ltla2num[fuseltla[lad21]]
        vax[0,d,i,a]+=dose1
        vax[1,d,i,a]+=dose2
  with open(fn,'wb') as fp:
    pickle.dump((vax,vaxday0),fp)

cumcases=np.cumsum(cases,0)

ncasedays=cases.shape[0];maxcaseday=caseday0+ncasedays
nvaxdays=vax.shape[1];maxvaxday=vaxday0+nvaxdays
vaxslice=np.mean(vax[:,maxcaseday-14-vaxeffecttime-vaxday0,:,:],0)# Average 1st+2nd doses
#vaxslice=vax[1,maxcaseday-14-vaxeffecttime-vaxday0,:,:]
cumcasesslice=cumcases[-14,:,:]#-cumcases[-100,:,:]
back=-0
week0=cases[-190-14:-190-7,:,:].sum(axis=0)
week1=cases[back-14:ncasedays+back-7,:,:].sum(axis=0)
week2=cases[back-7:ncasedays+back-0,:,:].sum(axis=0)

def err(ieff,veff,week1,week2,cumcasesslice,vaxslice,ltlapop):
  pred0=week1*(1-ieff/car*cumcasesslice/ltlapop)*(1-veff*vaxslice/ltlapop)
  G=(pred0*week2).sum()/(pred0*pred0).sum()
  #G=week2.sum()/pred0.sum()
  #G=1.4
  pred1=G*pred0
  e=pred1-week2
  return (e*e).sum()

if 0:
  a=0
  week1=week1.sum(axis=a)
  week2=week2.sum(axis=a)
  cumcasesslice=cumcasesslice.sum(axis=a)
  vaxslice=vaxslice.sum(axis=a)
  ltlapop=ltlapop.sum(axis=a)

ieff=1

for a in range(nages):
  print(num2age[a])
  for veff in [.1*i for i in range(11)]:
    print("%5.3f %5.3f %12g"%(ieff,veff,err(ieff,veff,week1[:,a],week2[:,a],cumcasesslice[:,a],vaxslice[:,a],ltlapop[:,a])))
  print()
  
if 0:
  a=8
  rat=(week2+1)/(week1+1)
  with open('tempg','w') as fp:
    for i in range(nltlas):
      print(cumcasesslice[i,a]/car/ltlapop[i,a],vaxslice[i,a]/ltlapop[i,a],rat[i,a],file=fp)
  print("Written graph file 'tempg' for age group",num2age[a])

if 0:
  #yy=np.log(week2.sum(axis=1)/week1.sum(axis=1))
  yy=np.log(week2.sum(axis=1))
  with open('tempg','w') as fp:
    for i in range(nltlas):
      print(sum(cumcasesslice[i,:])/car/sum(ltlapop[i,:]),sum(vaxslice[i,:])/sum(ltlapop[i,:]),yy[i],file=fp)
  print("Written graph file 'tempg'")

if 0:
  yy=np.log(week2)
  with open('tempg','w') as fp:
    for i in range(nltlas):
      for a in range(nages):
        print(cumcasesslice[i,a]/car/ltlapop[i,a],vaxslice[i,a]/ltlapop[i,a],yy[i,a],file=fp)
  print("Written graph file 'tempg'")

if 0:
  yy0=np.log(week0.sum(axis=1))
  yy1=np.log(week2.sum(axis=1))
  with open('tempg','w') as fp:
    for i in range(nltlas):
      print(sum(vaxslice[i,:])/sum(ltlapop[i,:]),yy0[i],yy1[i],file=fp)
  print("Written graph file 'tempg'")

if 1:
  a=8
  yy0=np.log(week0[:,a])
  yy1=np.log(week1[:,a])
  with open('tempg','w') as fp:
    for i in range(nltlas):
      print(sum(vaxslice[i,:])/sum(ltlapop[i,:]),yy0[i],yy1[i],file=fp)
  print("Written graph file 'tempg'")
