from stuff import *
import sys
from scipy.optimize import minimize
from math import log,exp,sqrt
import numpy as np
import re

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
def LondonLTLAs(ltla):
  "London LTLAs"
  return ltla2region[ltla]=='London'
def allLTLAs(x):
  "All LTLAs"
  return True
def BoltonLTLA(ltla):
  "Bolton LTLA"
  return ltla=='E08000001'

### Model ###
#
# Known:
# n_i      = number of confirmed cases on day i by specimen date (slightly adjusted for weekday)
# p        = case ascertainment rate (chance of seeing a case)
# r_j, s_j = Variant counts of non-B.1.617.2, B.1.617.2 in j^th week
# I_j      = set of days (week) corresponding to VOC counts r_j, s_j
# Assume chance of sequencing a case is a totally free parameter, and optimise over it
#
# Unknown:
# h   = daily growth advantage of B.1.617.2 over other variants
# g_i = growth rate of other variants on day i -> i+1
# A_0 = initial count of non-B.1.617.2
# B_0 = initial count of B.1.617.2
#
# Likelihood:
# A_{i+1}=e^{g_i}A_i
# B_{i+1}=e^{g_i+h}B_i
# n_i ~ Po(p(A_i+B_i))
# r_j ~ Po(q_j.A_{I_j})  (A_{I_j} means sum_{i in A_j}A_i)
# s_j ~ Po(q_j.B_{I_j})
# g_{i+1} ~ N(g_i,sig^2)
# q_j is optimised out, so q_j=(r_j+s_j)/(A_{I_j}+B_{I_j}), or effectively A_{I_j}/(A_{I_j}+B_{I_j}) ~ Beta(r_j+1,s_j+1)
# g_0 ~ N(g_{-1},v), where g_{-1} is from a case-based empirical calculation from a recent pre-B.1.617.2 period.
# [h ~ N(0,tau^2)]
#
### End Model ###


### Options ###

source="Sanger"
#source="COG-UK"
#source="SGTF"

# Can choose location size from "LTLA", "region", "country", "UK"
# Sanger works with LTLA, region, country
# COG-UK works with country, UK
# SGTF works with region, country
locationsize="LTLA"

exclude=set()
# exclude=set(['E08000001'])# This would exclude Bolton
include=allLTLAs
#include=LondonLTLAs
#include=BoltonLTLA

mgt=5# Mean generation time in days

# Earliest day to use case data
minday=datetoday('2021-04-01')# Inclusive

# Earliest day to use VOC count data, given as end-of-week. Will be rounded up to match same day of week as lastweek.
#firstweek=minday+6
firstweek=datetoday('2021-04-24')

nif1=0.5   # Non-independence factor for cases (less than 1 means downweight this information)
nif2=0.5   # Non-independence factor for VOC counts (ditto)
isd=1/0.05 # Inverse sd for prior on transmission advantage (as growth rate per day). 0 means uniform prior. 1/0.05 is fairly weak.
isd=1e-6   # Flat prior

# Effectively the prior on how much daily growth rate is allowed to change in 1 day
sig=0.002

# Case ascertainment rate
asc=0.4

# Discard this many cases at the end of the list of cases by specimen day
discarddays=2

# Collect together all locations without positive entries into one combined "Other" location
# (Makes little difference in practice)
bundleremainder=True

minopts={"maxiter":1000,"eps":1e-5}

### End options ###

print("Options")
print("Source:",source)
print("Location size:",locationsize)
print("Include:",include.__doc__)
print("Exclude:",exclude)
print("Generation time:",mgt,"days")
print("Earliest day for case data:",daytodate(minday))
print("Earliest week (using end of week date) to use VOC count data:",daytodate(firstweek))
print("nif1:",nif1)
print("nif2:",nif2)
print("Inverse sd for prior on growth:",isd)
print("Sigma (prior on daily growth rate change):",sig)
print("Case ascertainment rate:",asc)
print("Number of days of case data to discard:",discarddays)
print("Bundle remainder:",bundleremainder)
print("Minimiser options:",minopts)
print()

np.set_printoptions(precision=3,linewidth=120)

maxday=datetoday(max(apicases['date']))-discarddays# Inclusive
ndays=maxday-minday+1

if source=="Sanger":
  sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
  
  lastweek=datetoday(max(sanger['WeekEndDate']));assert maxday>=lastweek
  nweeks=(lastweek-firstweek)//7+1
  # Sanger week number is nweeks-1-(lastweek-day)//7

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
    if ltla in exclude or not include(ltla): continue
    day=datetoday(date)
    week=nweeks-1-(lastweek-day)//7
    if week>=0 and week<nweeks:
      place=reduceltla[ltla]
      if place not in vocnum: vocnum[place]=np.zeros([nweeks,2],dtype=int)
      if var=="B.1.617.2": vocnum[place][week][1]+=n
      else: vocnum[place][week][0]+=n
elif source=="COG-UK":
  cog=loadcsv("cog_metadata.csv")
  censor=6
  lastweek=datetoday(max(cog['sample_date']))-censor;assert maxday>=lastweek
  #nweeks=(lastweek-firstweek+1)//7
  nweeks=(lastweek-firstweek)//7+1
  # Week number is nweeks-1-(lastweek-day)//7
  
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
    week=nweeks-1-(lastweek-day)//7
    if week>=0 and week<nweeks:
      r=re.match("[^0-9-]*[0-9-]",seqname)
      coglab=seqname[:r.end()-1]
      place=reducecog(coglab)
      if place not in vocnum: vocnum[place]=np.zeros([nweeks,2],dtype=int)
      if var=="B.1.617.2": vocnum[place][week][1]+=1
      else: vocnum[place][week][0]+=1
elif source=="SGTF":
  sgtf=loadcsv("TechBriefing12Fig17.csv")
  lastweek=max(datetoday(x) for x in sgtf['week'])+6# Convert w/c to w/e convention
  assert maxday>=lastweek
  nweeks=(lastweek-firstweek)//7+1
  # Week number is nweeks-1-(lastweek-day)//7

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
    week=nweeks-1-(lastweek-day)//7
    if week>=0 and week<nweeks:
      place=reducesgtf(region)
      if place not in vocnum: vocnum[place]=np.zeros([nweeks,2],dtype=int)
      vocnum[place][week][int("SGTF" not in var)]+=n
  # Adjust for non-B.1.617.2 S gene positives, based on the assumption that these are in a non-location-dependent proportion to the number of B.1.1.7
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
      day=lastweek-7*(nweeks-1-week)
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
  if ltla not in reduceltla or ltla in exclude or not include(ltla): continue
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
#places.sort(key=lambda x: -vocnum[x].sum())# Descending order of voc counts
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
        x=(hmin+(i+.5)/ndiv*(hmax-hmin))*7# Convert to weekly growth rate
        a,b,c,d=M[0,0],M[0,1],M[1,0],M[1,1]
        l0=d*x-(betaln(a,b)+betaln(c,d))
        e=exp(x)
        res=quad(lambda z: exp( (b+d-1)*log(z) - (a+b)*log(1+z) - (c+d)*log(1+e*z) + l0 ), 0, inf)
        logp[i]+=log(res[0])
  g=log(tot[0,0]*tot[1,1]/(tot[0,1]*tot[1,0]))/7
  dg=sqrt((1/tot.flatten()).sum())/7
  print("Overall cross ratio:",Rdesc(g,dg),tot.flatten())
  print("Inverse variance weighting method using log(CR):",Rdesc(L0/L1/7,1/sqrt(L1)/7))
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
  
for w in range(nweeks-1):
  day0=lastweek-(nweeks-w)*7+1
  print(daytodate(day0),"-",daytodate(day0+13))
  crossratiosubdivide(vocnum[place][w:w+2] for place in places)
print("All week pairs")
crossratiosubdivide(vocnum[place][w:w+2] for place in places for w in range(nweeks-1))

print("Estimating transmission advantage using variant counts together with case counts")
print("================================================================================")
print()

# ndays+2 parameters to be optimised:
# 0: a0
# 1: b0
# 2: h/sig
# 3 ... 3+ndays-2 : g_0/sig, ..., g_{ndays-2}/sig
# (would be tidier to put sig in the likelihood instead of on the variables, but SLSQP seems
#  to get in trouble if its variables to be optimised are on too small a scale)

def expand(xx,sig):
  (a0,b0,h)=xx[:3]
  AA=[exp(a0)];BB=[exp(b0)]
  a=a0;b=b0
  for i in range(ndays-1):
    g=xx[3+i]
    a+=g*sig;b+=(g+h)*sig
    AA.append(exp(a))
    BB.append(exp(b))
  return AA,BB

# Return negative log likelihood
def NLL(xx,lcases,lvocnum,sig,p,lprecases):
  a,b=lprecases[0]+.5,lprecases[1]+.5
  g0=log(b/a)/7/sig
  v0=(1/a+1/b)/49/sig**2+1
  tot=-(xx[3]-g0)**2/(2*v0)
  AA,BB=expand(xx,sig)
  # Component of likelihood due to number of confirmed cases seen
  for i in range(ndays):
    lam=p*(AA[i]+BB[i])
    # max with -10000 because the expression is unbounded below which can cause a problem for SLSQP
    tot+=max((lcases[i]-lam+lcases[i]*log(lam))*nif1,-10000)
  # Term to regulate change in growth rate
  for i in range(ndays-2):
    # Could downweight (allow larger) changes in growth on or near roadmap days, but in practice that makes almost no difference
    tot+=-(xx[3+i+1]-xx[3+i])**2/2
  # Term to align the variant numbers with VOC count data
  for w in range(nweeks):
    endweek=lastweek-(nweeks-1-w)*7-minday
    A=sum(AA[endweek-6:endweek+1])
    B=sum(BB[endweek-6:endweek+1])
    tot+=(lvocnum[w][0]*log(A/(A+B))+lvocnum[w][1]*log(B/(A+B)))*nif2
  # Prior on h
  tot+=-(xx[2]*sig*isd)**2/2
  return -tot

def getlikelihoods(fixedh=None):
  summary={}
  logp=np.zeros(ndiv)
  TAA=np.zeros(ndays)
  TBB=np.zeros(ndays)
  for place in places:
    print(place)
    print("="*len(place))
    print()
    print("                        Nonvar    Var   Seen")
    for w in range(nweeks):
      day0,day1=lastweek-(nweeks-w)*7+1,lastweek-(nweeks-1-w)*7
      print(daytodate(day0),"-",daytodate(day1),"%6d %6d %6.0f"%(vocnum[place][w][0],vocnum[place][w][1],sum(cases[place][day0-minday:day1-minday+1])))
    print()
    xx=np.zeros(ndays+2)
    bounds=[(-10,20),(-10,20),(-1/sig,1/sig)]+[(-1/sig,1/sig)]*(ndays-1)
    if fixedh!=None: xx[2]=fixedh;bounds[2]=(fixedh,fixedh)
    res=minimize(NLL,xx,args=(cases[place],vocnum[place],sig,asc,precases[prereduce(place)]),bounds=bounds,method="SLSQP",options=minopts)
    if not res.success: raise RuntimeError(res.message)
    xx=res.x
    AA,BB=expand(xx,sig)
    TAA+=AA;TBB+=BB
    (a0,b0,h)=xx[:3]
    print("A    = estimated number of new cases of non-B.1.617.2 on this day multiplied by the ascertainment rate")
    print("B    = estimated number of new cases of B.1.617.2 on this day multiplied by the ascertainment rate")
    print("Pred = predicted number of cases seen this day = A+B")
    print("Seen = number of cases seen this day, after weekday adjustment")
    print("Q    = estimated reproduction rate of non-B.1.617.2 on this day")
    print("R    = estimated reproduction rate of B.1.617.2 on this day")
    print()
    print("      Date       A       B    Pred    Seen       Q     R")
    for i in range(ndays):
      print(daytodate(minday+i),"%7.0f %7.0f %7.0f %7.0f"%(asc*AA[i],asc*BB[i],asc*(AA[i]+BB[i]),cases[place][i]),end='')
      if i<ndays-1:
        g=xx[3+i]
        Q,R=(exp(g*sig*mgt),exp((g+h)*sig*mgt))
        print("   %5.2f %5.2f"%(Q,R))
      else:
        print()
    print()
    if fixedh!=None: summary[place]=(Q,R);continue
    h0=xx[2]
    ff=[0,res.fun,0]
    eps=0.01/sig
    for i in [-1,1]:
      h=h0+i*eps
      xx=[0,0,h]+[0]*(ndays-1)
      bounds[2]=(h,h)
      res=minimize(NLL,xx,args=(cases[place],vocnum[place],sig,asc,precases[prereduce(place)]),bounds=bounds,method="SLSQP",options=minopts)
      if not res.success:
        print(res)
        print(place)
        print("xx =",xx)
        print("lcases =",list(cases[place]))
        print("lprecases =",precases[prereduce(place)])
        print("lvocnum =",vocnum[place])
        print("sig =",sig)
        print("asc =",asc)
        print("bounds =",bounds)
        print("nif1 =",nif1)
        print("nif2 =",nif2)
        print("nweeks, ndays, minday, lastweek =",nweeks,",",ndays,",",minday,",",lastweek)
        print("isd =",isd)
        print("minopts =",minopts)
        raise RuntimeError(res.message)
      ff[i+1]=res.fun
    # Use observed Fisher information to make confidence interval
    fi=(ff[0]-2*ff[1]+ff[2])/eps**2
    if fi>0:
      dh=1.96/sqrt(fi)
    else:
      dh=100/sig
    (Tmin,T,Tmax)=[(exp(h*sig*mgt)-1)*100 for h in [h0-dh,h0,h0+dh]]
    print("Estimated transmission advantage %.0f%% (%.0f%% - %.0f%%)"%(T,Tmin,Tmax))
    summary[place]=(Q,R,Tmin,T,Tmax)
    print()
    print("    g     T    log lik")
    for i in range(ndiv):
      h=(hmin+(hmax-hmin)*i/(ndiv-1))/sig
      xx=[0,0,h]+[0]*(ndays-1)
      bounds[2]=(h,h)
      res=minimize(NLL,xx,args=(cases[place],vocnum[place],sig,asc,precases[prereduce(place)]),bounds=bounds,method="SLSQP",options=minopts)
      if not res.success: raise RuntimeError(res.message)
      logp[i]+=ff[1]-res.fun
      print("%5.3f %5.3f  %9.2f"%(h*sig,exp(h*sig*mgt),logp[i]))
    print()
    sys.stdout.flush()
  print()
  return summary,logp,TAA,TBB

ndiv=11
hmin=0.03;hmax=0.15

summary,logp,TAA,TBB=getlikelihoods()

print("Location                       Q     R      T")
for place in places:
  (Q,R,Tmin,T,Tmax)=summary[place]
  print("%-25s  %5.2f %5.2f  %4.0f%% ( %4.0f%% - %4.0f%% )"%(place,Q,R,T,Tmin,Tmax))
print()
print("Q = point estimate of reproduction rate of non-B.1.617.2 on",daytodate(maxday-1))
print("R = point estimate of reproduction rate of B.1.617.2 on",daytodate(maxday-1))
print("T = estimated transmission advantage = R/Q as a percentage increase")
print()

print("    g     T    log lik")
for i in range(ndiv):
  g=(hmin+(hmax-hmin)*i/(ndiv-1))
  print("%5.3f %5.3f  %9.2f"%(g,exp(g*mgt),logp[i]))
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
print("Combined growth advantage per day %.3f (%.3f - %.3f)"%(h0,h0-dh,h0+dh))
print("Combined transmission advantage %.0f%% (%.0f%% - %.0f%%) (assuming fixed generation time of %g days)"%(T,Tmin,Tmax,mgt))
print()

print("Re-running using global optimum growth advantage")
print()
summary,logp,TAA,TBB=getlikelihoods(fixedh=h0/sig)

print("Total predicted counts using global optimum growth advantage")
print()

print("A    = estimated number of new cases of non-B.1.617.2 on this day multiplied by the ascertainment rate")
print("B    = estimated number of new cases of B.1.617.2 on this day multiplied by the ascertainment rate")
print("Pred = predicted number of cases seen this day = A+B")
print("Seen = number of cases seen this day, after weekday adjustment")
print("Q    = estimated reproduction rate of non-B.1.617.2 on this day")
print("R    = estimated reproduction rate of B.1.617.2 on this day")
print()
print("      Date       A       B    Pred    Seen       Q     R")
for i in range(ndays):
  print(daytodate(minday+i),"%7.0f %7.0f %7.0f %7.0f"%(asc*TAA[i],asc*TBB[i],asc*(TAA[i]+TBB[i]),sum(cases[place][i] for place in places)),end='')
  if i<ndays-1:
    g=log(TAA[i+1]/TAA[i])
    Q,R=(exp(g*mgt),exp((g+h0)*mgt))
    print("   %5.2f %5.2f"%(Q,R))
  else:
    print()
print()

print("Location                       Q     R")
for place in places:
  (Q,R)=summary[place]
  print("%-25s  %5.2f %5.2f"%(place,Q,R))
print()
print("Q = point estimate of reproduction rate of non-B.1.617.2 on",daytodate(maxday-1))
print("R = point estimate of reproduction rate of B.1.617.2 on",daytodate(maxday-1))
print()
print("Combined growth advantage per day %.3f (%.3f - %.3f)"%(h0,h0-dh,h0+dh))
print("Combined transmission advantage %.0f%% (%.0f%% - %.0f%%) (assuming fixed generation time of %g days)"%(T,Tmin,Tmax,mgt))
