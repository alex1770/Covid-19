from stuff import *
from math import log,exp
from scipy.optimize import minimize
from random import randrange
from math import sqrt
import sys
import numpy as np
np.set_printoptions(precision=6,linewidth=250,suppress=True)

minday=Date('2021-11-25')
maxday=Date('2021-12-25')# Only go up to dates strictly before this one
pubday=getpublishdate()
discarddays=3# Discard last few case counts by specimen date since these are incomplete (irrelevant here because we're stopping much earlier anyway)
mincount=5
step=7
outdir='gentimeoutput'
os.makedirs(outdir,exist_ok=True)
conf=0.95
adjustbycases=True
print("Adjust by cases:",adjustbycases)
nsamp=1000
if len(sys.argv)>1: nsamp=int(sys.argv[1])
mode="byage"
print("Mode",mode)

ages=[(a,a+10) for a in range(0,80,10)]+[(80,150)]
nages=len(ages)

ONSpop={}
for (desc,acode,sex,age,n) in csvrows('ONS-population_2021-08-05.csv',['category','areaCode','gender','age','population']):
  if desc=='AGE_SEX_5YEAR' and sex=='ALL' and acode=='E92000001':# England
    pa=parseage(age)
    for (a0,a1) in ages:
      if pa[0]>=a0 and pa[1]<=a1: ONSpop[(a0,a1)]=ONSpop.get((a0,a1),0)+int(n)

sp0,sp=getcasesbyagespeccomplete(minday=minday,maxday=pubday,ages=ages,location='England')
casesbyregion={ages[a]:sp[:pubday-discarddays-minday,a] for a in range(nages)}
casesbyregion['England']=sum(casesbyregion.values())
nspec=casesbyregion['England'].shape[0]

# From fig. 8B of Tech Briefing 33, https://www.gov.uk/government/publications/investigation-of-sars-cov-2-variants-technical-briefings
avd=loadcsv('age_variant_series.csv')
navd=max(Date(x+'US') for x in avd['spec.date'])+1-minday
vocnum={age:np.zeros([navd,2],dtype=int) for age in ages}
for (date,var,age,n) in zip(avd['spec.date'],avd['variant'],avd['age.group'],avd['count_age']):
  d=Date(date+'US')-minday
  if d>=0:
    a=parseage(age)
    l=[b for b in ages if a[0]>=b[0] and a[1]<=b[1]]
    assert len(l)<=1
    if len(l)==1:
      vocnum[l[0]][d,int(var=='Omicron')]+=n
nsgtf=navd
vocnum['Allages']=sum(vocnum.values())

n=min(nsgtf,nspec)
cv={}
for loc in vocnum:
  if loc=='Allages': continue
  vv=vocnum[loc][:n,:]
  if adjustbycases:
    cases=casesbyregion[loc][:n]
    sprop=vv/(vv.sum(axis=1)[:,None])
    cv[loc]=cases[:,None]*sprop
  else:
    cv[loc]=vv

def muggle(age):
  if age[1]<150: return "%d-%d"%(age[0],age[1]-1)
  return "%d+"%age[0]
    
with open('age-variant.csv','w') as fp:
  print("Date,",end='',file=fp)
  for age in ages:
    m=muggle(age)
    print("Delta_%s,Omicron_%s"%(m,m),end='',file=fp)
    if age!=ages[-1]: print(",",end='',file=fp)
  print(file=fp)
  for d in range(n):
    print(minday+d+',',end='',file=fp)
    for age in ages:
      print("%.3f,%.3f"%(cv[age][d][0]/ONSpop[age]*1e5,cv[age][d][1]/ONSpop[age]*1e5),end='',file=fp)
      if age!=ages[-1]: print(',',end='',file=fp)
    print(file=fp)

