from stuff import *
import sys,re,argparse,pickle,csv
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.special import gammaln
from math import log,exp,sqrt,sin,pi
import numpy as np
from subprocess import Popen,PIPE
from datetime import datetime

# ltla_age.csv from https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=newCasesBySpecimenDateAgeDemographics&format=csv
# wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=vaccinationsAgeDemographics&format=csv' -O vaccinedata.csv
# Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv

# NCU:
# COG-UK data from https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv
# SGTF   data from Fig.16 Tech Briefing 12: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/988608/Variants_of_Concern_Technical_Briefing_12_Data_England.xlsx

# Case age bands (all English LTLAs): 0-5, 5-10, ..., 85-90, 90+ (18 5-year bands and 90+)
# Vax and NIMS-population age bands for English LTLAs: 18-25, 25-30, 30-35, ..., 85-90, 90+ (18-25, 13 5-year bands, 90+)
# Standardise on these age bands: 0-5, 5-10, 10-15, 15-18, 18-25, 25-30, ..., 85-90, 90+
# by (i) zero-extending NIMS data, and
#    (ii) using simple interpolation to convert case data from 15-20, 20-25 to 15-18, 18-25
# These bands are hard-coded, perforce, including the assumption that there are equal number of case age bands and unified age bands (19)
num2caseage=['%02d_%02d'%(i,i+4) for i in range(0,90,5)]+['90+']
caseage2num=dict((a,i) for (i,a) in enumerate(num2caseage))
caseage2num['00_59']=None
caseage2num['60+']=None
caseage2num['unassigned']=None

# Unified age bands; serves for decoding vax data too.
num2age=['00_05','05_10','10_15','15_18','18_24']+['%02d_%02d'%(i,i+4) for i in range(25,90,5)]+['90+']
age2num=dict((a,i) for (i,a) in enumerate(num2age))
numage=len(num2age)

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
# 2. Dashboard/api use LAD19 except that E09000001 (City of London) has been fused into E09000012 (Hackney).
#    and E06000053 (Isles of Scilly) is fused into E06000052 (Cornwall).
# 3. NIMS (population and vaccine data) prior to July 2021 uses LAD20 (cf LAD19: E0700000[4567] fused into E06000060).
# 4. NIMS (population and vaccine data) from July 2021 uses LAD21 (cf LAD20: E0700015[0236] fused into E06000061, E0700015[145] fused into E06000062).
# Easier to fuse than to split, so standardise on LAD21 with E09000001 fused into E09000012 and E06000053 fused into E06000052.
fuseltla=dict(zip(ltlaengdata['LAD21CD'],ltlaengdata['LAD21CD']))
fuseltla['E07000004']='E06000060'
fuseltla['E07000005']='E06000060'
fuseltla['E07000006']='E06000060'
fuseltla['E07000007']='E06000060'
fuseltla['E09000001']='E09000012'
fuseltla['E06000053']='E06000052'
fuseltla['E07000150']='E06000061'
fuseltla['E07000152']='E06000061'
fuseltla['E07000153']='E06000061'
fuseltla['E07000156']='E06000061'
fuseltla['E07000151']='E06000062'
fuseltla['E07000154']='E06000062'
fuseltla['E07000155']='E06000062'
# Convert standardised (fused) ltlas to and from numbers:
num2ltla=sorted(list(set(fuseltla.values())))
ltla2num=dict((ltla,i) for (i,ltla) in enumerate(num2ltla))
numltla=len(num2ltla)
s=set(ltlaengdata['LAD21CD']);assert all(x in s for x in fuseltla.values())

indir='inputs'
os.makedirs(indir,exist_ok=True)
now=datetime.utcnow()
apiday=datetoday(now.strftime('%Y-%m-%d'))-(now.hour<16)

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
asc=0.4

# Discard this many cases at the end of the list of cases by specimen day
discardcasedays=2# Make need to increase to 3 to allow for Wales and NI late reporting at weekends and bank holidays
missingcasedays=4# Expect 4 missing days when using newCasesBySpecimenDateAgeDemographics; setting missingcasedays at 4 means it will refresh cases every day
missingvaxdays=12# Only expect 0 missing days using vaccinationsAgeDemographics, but we only need to refresh every 12 days or so because vaccines don't have an immediate effect

minopts={"maxiter":10000,"eps":1e-4}

voclen=(1 if source=="COG-UK" else 7)

conf=0.95

### End options ###

opts={
  "Source": source,
  "LTLA set": ltlaset,
  "LTLA exclude": list(ltlaexclude),
  "Generation time (days)": mgt,
  "Earliest week (using end of week date) to use VOC count data": daytodate(firstweek),
  "Case ascertainment rate": asc,
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

source=opts["Source"]
ltlaset=opts["LTLA set"]
ltlaexclude=set(opts["LTLA exclude"])
mgt=opts["Generation time (days)"]
firstweek=datetoday(opts["Earliest week (using end of week date) to use VOC count data"])
asc=opts["Case ascertainment rate"]
discardcasedays=opts["Number of days of case data to discard"]
minopts=opts["Minimiser options"]
voclen=opts["Length of time period over which VOC counts are given (days)"]
conf=opts["Confidence level"]

zconf=norm.ppf((1+conf)/2)

np.set_printoptions(precision=3,linewidth=150)

print("Options:")
print()
for x in sorted(list(opts)): print("%s:"%x,opts[x])
print()
sys.stdout.flush()

okplaces=set(ltla for ltla in fuseltla.values() if includeltla(ltla,ltlaset) and not ltla in ltlaexclude)
places=sorted(list(okplaces))

npl=len(places)
# ltla2num={ltla:i for (i,ltla) in enumerate(places)}

regions=sorted(list(set(ltla2region.values())))
nreg=len(regions)
reg2num={reg:i for (i,reg) in enumerate(regions)}
ltla2regnum=np.array([reg2num[ltla2region[ltla]] for ltla in places])

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

# Alter
ltlapopdata=loadcsv("LTLA-NIMS-populations.LAD21.18+.csv")
ltlapop=np.zeros([numltla,numage],dtype=int)
regionpop=np.zeros(nreg,dtype=int)
for (lad21,n1,n2) in zip(ltlapopdata['LTLA Code'],ltlapopdata['Under 16'],ltlapopdata['16+']):
  ltla=fuseltla[lad21]
  if ltla not in okplaces: continue
  i=ltla2num[ltla]
  ltlapop[i]+=n1+n2
  r=ltla2regnum[i]
  regionpop[r]+=n1+n2
rpopratio=np.array([ltlapop[i]/regionpop[ltla2regnum[i]] for i in range(numltla)])
tpopratio=ltlapop/ltlapop.sum()

# Input is of the form "5_9", "25_29", "60+" etc
def parseage_api(age):
  if '_' in age: x=age.split('_');return int(x[0]),int(x[1])+1
  if age[-1]=='+': return int(age[:-1]),150
  if age=='unassigned': return 0,150
  else: raise RuntimeError("Unrecognised age band: "+age)

# Input is of the form: D*.0-50, D*.50-55, ..., D1.0-30, D1.30-35, D1.35-40, D1.40-45, ..., D2.40-45, D2.45-50, ...
def parseage_nimsvax(age):
  if age[:1]!='D': return None,None
  f=age.find('.')
  g=age.find('-',f+1)
  return (age[:f],(int(age[f+1:g]),int(age[g+1:])))

# Under 16,16-29,30-34,35-39,40-44,45-49,50-54,55-59,60-64,65-69,70-74,75-79,80+,16+
def parseage_nimspop(age):
  if '-' in age: x=age.split('-');return int(x[0]),int(x[1])+1
  if age[:5]=='Under': return 0,int(age[6:])
  if age[-1]=='+': return int(age[:-1]),150
  return None

# Get case data into caseages[daynum,ltlanum,agebandnum]
# The saved (pickled) format of case ages uses standardised (fused) ltlas, and standardised age bands
fn=os.path.join(indir,'ltla_age_cases.pickle')
ok=0
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    caseages,caseday0=pickle.load(fp)
  ndays,nltlas,nages=caseages.shape
  if caseday0+ndays>=apiday-missingcasedays: ok=1
if not ok:
  print("Loading case-age data from api")
  la=loadcsv_it(api_v2('areaType=ltla&metric=newCasesBySpecimenDateAgeDemographics&format=csv').text.split('\n'))
  print("Processing case-age data")
  caseday0=datetoday(min(la['date']))
  caseday1=datetoday(max(la['date']))
  n=caseday1-caseday0+1
  caseages=np.zeros([n,numltla,len(num2caseage)])
  for lad19,date,age,cases in zip(la['areaCode'],la['date'],la['age'],la['cases']):
    a=caseage2num[age]
    if a!=None:
      ltla=fuseltla[lad19]
      caseages[datetoday(date)-caseday0,ltla2num[ltla],caseage2num[age]]=cases
  t=caseages[:,:,3]+caseages[:,:,4]# (15-20) + (20-25)
  caseages[:,:,3]=3/10*t# Estimate of 15-18
  caseages[:,:,4]=7/10*t# Estimate of 18-25
  with open(fn,'wb') as fp:
    pickle.dump((caseages,caseday0),fp)

# The saved (pickled) format of case ages uses standardised (fused) ltlas, and standardised age bands
fn=os.path.join(indir,'vax_age_cases.pickle')
ok=0
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    vaxages,vaxday0=pickle.load(fp)
  ndays,nltlas,nages=vaxages.shape
  if vaxday0+ndays>=apiday-missingvaxdays: ok=1
if not ok:
  print("Loading vax-age data from api")
  #la=loadcsv_it(api_v2('areaType=ltla&metric=vaccinationsAgeDemographics&format=csv').text.split('\n'))
  la=loadcsv("vaxdata.csv")
  print("Processing vax-age data")
  vaxday0=datetoday(min(la['date']))
  vaxday1=datetoday(max(la['date']))
  n=vaxday1-vaxday0+1
  vax1=np.zeros([n,numltla,numage])
  vax2=np.zeros([n,numltla,numage])
  for lad21,date,age,dose1,dose2 in zip(la['areaCode'],la['date'],la['age'],la['cumPeopleVaccinatedFirstDoseByVaccinationDate'],la['cumPeopleVaccinatedSecondDoseByVaccinationDate']):
    a=age2num[age]
    if a!=None:
      ltla=fuseltla[lad21]
      vaxages[datetoday(date)-vaxday0,ltla2num[ltla],vaxage2num[age]]=vaxs
  t=vaxages[:,:,3]+vaxages[:,:,4]# (15-20) + (20-25)
  vaxages[:,:,3]=3/10*t# Estimate of 15-18
  vaxages[:,:,4]=7/10*t# Estimate of 18-25
  with open(fn,'wb') as fp:
    pickle.dump((vaxages,vaxday0),fp)


abcd

from random import random,seed
#seed(42)
pvax={place:[[] for w in range(nweeks-1)] for place in places}
for w in range(nweeks-1):
  date=daytodate(firstweek+w*7+10-vaxeffecttime)
  id=max(dt for dt in vaxdat if dt<date)
  print("Using NIMS vax w/e %s to correspond to Sanger w/e %s -> w/e %s"%(id,daytodate(firstweek+w*7),daytodate(firstweek+w*7+7)))
  v=vaxdat[id]
  
  vaxnum={}# Map from LTLA -> age -> numvaxed
  for (i,lad21) in enumerate(v['LTLA Code']):
    ltla=fuseltla[lad21]
    if ltla in okplaces:
      if ltla not in vaxnum: vaxnum[ltla]={}
      for age in v:
        d,a=parseage_nimsvax(age)
        if a!=None:
          vaxnum[ltla][a]=vaxnum[ltla].get(a,0)+v[age][i]

  vaxpop={}# Map from LTLA -> age -> population
  for (i,lad21) in enumerate(ltlapopdata['LTLA Code']):
    ltla=fuseltla[lad21]
    if ltla in okplaces:
      if ltla not in vaxpop: vaxpop[ltla]={}
      for age in ltlapopdata:
        a=parseage_nimspop(age)
        if a!=None:
          vaxpop[ltla][a]=vaxpop[ltla].get(a,0)+ltlapopdata[age][i]
  
  for ltla in caseages:
    if ltla in okplaces:
      cas=caseages[ltla]
      vax=vaxnum[ltla]
      pop=vaxpop[ltla]
      pp=[]
      for ca in cas:
        n=cas[ca][w+1]
        if n==0: continue
        num=den=0
        for va in vax:
          # Want (vax number) as weighted by P(case age interval|vax age interval) = |ca intersect va|/|va|
          num+=max(min(va[1],ca[1],100)-max(va[0],ca[0]),0)/(min(va[1],100)-va[0])*vax[va]
        for pa in pop:
          den+=max(min(pa[1],ca[1],100)-max(pa[0],ca[0]),0)/(min(pa[1],100)-pa[0])*pop[pa]
        pp.append((n,min(num/den,1)))
      pvax[ltla][w]=pp

xx,L=optimise()
print()

print("Variables:",xx)
print("Log likelihood:",L)
NLL(xx*condition,const=True,pic=True)
H=Hessian(xx)
print("Hessian:");print(H)
eig=np.linalg.eigh(H)
print("Eigenvalues:",eig[0])
print()

h=xx[0]/7;dh=1/sqrt(H[0,0])/7
print("Logarithmic growth rate advantage/day: %.1f%% (%.1f%% - %.1f%%)"%(h*100,(h-zconf*dh)*100,(h+zconf*dh)*100))
print("Multiplicative growth rate advantage/day: %.1f%% (%.1f%% - %.1f%%)"%((exp(h)-1)*100,(exp(h-zconf*dh)-1)*100,(exp(h+zconf*dh)-1)*100))
print("R-number advantage: %.2f (%.2f - %.2f)"%(exp(mgt*h),exp(mgt*(h-zconf*dh)),exp(mgt*(h+zconf*dh))))
print()

print("Confidence intervals from single variables:")
for i in range(N):
  x=xx[i];dx=1/sqrt(H[i,i])
  print("%10s: %5.3f (%5.3f - %5.3f)"%(params[i][0],x,x-zconf*dx,x+zconf*dx))
print()

print("Confidence intervals from multivariate calculation:")
C=np.linalg.inv(H)
for i in range(N):
  x=xx[i];dx=sqrt(C[i,i])
  print("%10s: %5.3f (%5.3f - %5.3f)"%(params[i][0],x,x-zconf*dx,x+zconf*dx))
print()

samp,cc=makesamples(xx,H)
print("Confidence intervals from multivariate simulation:")
for i in range(N):
  print("%10s: %5.3f (%5.3f - %5.3f)"%((params[i][0],)+cc[i]))

if 0:
  l=list(pvax)
  l.sort(key=lambda x: pvax[x][3])
  for x in l:
    print(x,"%5.0f %5.0f %5.0f %5.0f"%tuple([z*100 for z in pvax[x]]),"  ",ltla2name[x])
