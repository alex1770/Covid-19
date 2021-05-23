from stuff import *
import sys
from scipy.optimize import minimize
from math import log,exp,sqrt
import numpy as np
import re

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
# h ~ N(0,tau^2)
#
### End Model ###


### Options ###

exclude=set()
# exclude=set(['E08000001'])# This would exclude Bolton

#source="Sanger"
#source="COG-UK"
source="SGTF"

mgt=5# Mean generation time in days

# Earliest day to use case data
minday=datetoday('2021-04-01')# Inclusive

# Earliest day to use VOC count data, given as end-of-week. Will be rounded up to match same day of week as lastweek.
#firstweek=minday+6
firstweek=datetoday('2021-04-24')

# Can only use these options with Sanger or SGTF data
#reduceltla=ltla2ltla
reduceltla=ltla2region
reducesgtf=sgtf2region

#reduceltla=ltla2country
#reducecog=coglab2country
#reducesgtf=sgtf2country

#reduceltla=ltla2uk
#reducecog=coglab2uk

#reducecog=coglab2coglab# Can't use this very accurately because there isn't a consistent LTLA->coglab map

nif1=0.5   # Non-independence factor for cases (less than 1 means downweight this information)
nif2=0.5   # Non-independence factor for VOC counts (ditto)
isd=1/0.05 # Inverse sd for prior on transmission advantage (as growth rate per day). 0 means uniform prior. 1/0.05 is fairly weak.

# Effectively the prior on how much daily growth rate is allowed to change in 1 day
sig=0.002

# Case ascertainment rate
asc=0.4

discarddays=2

### End options ###

np.set_printoptions(precision=3,linewidth=120)

maxday=datetoday(max(apicases['date']))-discarddays# Inclusive
ndays=maxday-minday+1

if source=="Sanger":
  sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
  
  lastweek=datetoday(max(sanger['WeekEndDate']));assert maxday>=lastweek
  nweeks=(lastweek-firstweek)//7+1
  # Sanger week number is nweeks-1-(lastweek-day)//7
  
  # Get Sanger (variant) data into a suitable form
  vocnum={}
  for (date,ltla,var,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
    if ltla in exclude: continue
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
  # 2021-02-04  13258     765   5.77%
  # 2021-02-11  16540     681   4.12%
  # 2021-02-18  15707     537   3.42%
  # 2021-02-25  16676     472   2.83%
  # 2021-03-04  14524     372   2.56%
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
cases={}
for (ltla,date,n) in zip(apicases['areaCode'],apicases['date'],apicases['newCasesBySpecimenDate']):
  if ltla not in reduceltla or ltla in exclude: continue
  day=datetoday(date)
  d=day-minday
  if d<0 or d>=ndays: continue
  place=reduceltla[ltla]
  if place not in vocnum: continue
  if place not in cases: cases[place]=np.zeros(ndays)
  cases[place][d]+=n/weekadjp[day%7]
places=sorted(list(cases))
#for x in places: print(x);print(cases[x]);print()

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
def NLL(xx,lcases,lvocnum,sig,p):
  #print(xx)
  AA,BB=expand(xx,sig)
  tot=0
  # Component of likelihood due to number of confirmed cases seen
  for i in range(ndays):
    lam=p*(AA[i]+BB[i])
    tot+=(-lam+lcases[i]*log(lam))*nif1
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

summary={}
for place in places:
  print(place)
  print("="*len(place))
  print()
  print("                        Nonvar    Var   Seen")
  for w in range(nweeks):
    day0,day1=lastweek-(nweeks-w)*7+1,lastweek-(nweeks-1-w)*7
    print(daytodate(day0),"-",daytodate(day1),"%6d %6d %6.0f"%(vocnum[place][w][0],vocnum[place][w][1],sum(cases[place][day0-minday:day1-minday+1])))
  print()
  bounds=[(-10,20),(-10,20),(-1/sig,1/sig)]+[(-1/sig,1/sig)]*(ndays-1)
  res=minimize(NLL,[0]*(ndays+2),args=(cases[place],vocnum[place],sig,asc),bounds=bounds,method="SLSQP",options={"maxiter":1000})
  if not res.success: raise RuntimeError(res.message)
  xx=res.x
  #print(res.fun)
  #print(xx)
  AA,BB=expand(xx,sig)
  (a0,b0,h)=xx[:3]
  print("A    = estimated number of new cases of non-B.1.617.2 on this day")
  print("B    = estimated number of new cases of B.1.617.2 on this day")
  print("Pred = predicted number of cases seen this day = (ascertainment rate)*(A+B)")
  print("Seen = number of cases seen this day, after weekday adjustment")
  print("Q    = estimated reproduction rate of non-B.1.617.2 on this day")
  print("R    = estimated reproduction rate of B.1.617.2 on this day")
  print()
  print("      Date       A       B    Pred    Seen       Q     R")
  for i in range(ndays):
    print(daytodate(minday+i),"%7.0f %7.0f %7.0f %7.0f"%(AA[i],BB[i],asc*(AA[i]+BB[i]),cases[place][i]),end='')
    if i<ndays-1:
      g=xx[3+i]
      Q,R=(exp(g*sig*mgt),exp((g+h)*sig*mgt))
      print("   %5.2f %5.2f"%(Q,R))
    else:
      print()
  h0=xx[2]
  ff=[0,res.fun,0]
  eps=0.01/sig
  for i in [-1,1]:
    h=h0+i*eps
    xx=[0,0,h]+[0]*(ndays-1)
    bounds[2]=(h,h)
    res=minimize(NLL,xx,args=(cases[place],vocnum[place],sig,asc),bounds=bounds,method="SLSQP",options={"maxiter":1000})
    if not res.success: raise RuntimeError(res.message)
    ff[i+1]=res.fun
  # Use observed Fisher information to make confidence interval
  dh=1.96/sqrt((ff[0]-2*ff[1]+ff[2])/eps**2)
  (Tmin,T,Tmax)=[(exp(h*sig*mgt)-1)*100 for h in [h0-dh,h0,h0+dh]]
  print("Transmission advantage %.0f%% (%.0f%% - %.0f%%)"%(T,Tmin,Tmax))
  summary[place]=(Q,R,Tmin,T,Tmax)
  print()
print()

print("Location                       Q     R      T")
for place in places:
  (Q,R,Tmin,T,Tmax)=summary[place]
  print("%-25s  %5.2f %5.2f  %4.0f%% ( %4.0f%% - %4.0f%% )"%(place,Q,R,T,Tmin,Tmax))
print()
print("Q = point estimate of reproduction rate of non-B.1.617.2 on",daytodate(maxday-1))
print("R = point estimate of reproduction rate of B.1.617.2 on",daytodate(maxday-1))
print("T = estimated transmission advantage = R/Q as a percentage increase")
