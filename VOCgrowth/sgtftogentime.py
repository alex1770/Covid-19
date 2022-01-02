from stuff import *
from math import log,exp
from scipy.optimize import minimize
import numpy as np
np.set_printoptions(precision=3,linewidth=250,suppress=True)

l=[x for x in os.listdir('.') if x[:19]=='sgtf_regionepicurve']
if l==[]: raise RuntimeError("No sgtf_regionepicurve*.csv file found in current directory; download from https://www.gov.uk/government/publications/covid-19-omicron-daily-overview")
sgtf=loadcsv(max(l))
regions=sorted(list(set(sgtf['UKHSA_region'])))
minday=Date('2021-11-25')
maxday=Date('2021-12-25')# Only go up to dates strictly before this one
pubday=getpublishdate()
discard=3# Discard last few case counts by specimen date since these are incomplete (irrelevant here because we're stopping much earlier anyway)

ages=[(0,150)]
casesbyregion={}
for loc in regions:
  sp0,sp=getcasesbyagespeccomplete(minday=minday,maxday=pubday,ages=ages,location=loc)
  casesbyregion[loc]=sp[:pubday-discard-minday]
casesbyregion['England']=sum(casesbyregion.values())
nspec=casesbyregion['England'].shape[0]
  
# Get SGTF data into a suitable form
vocnum={}
background=[0,0]
nsgtf=max(Date(d) for d in sgtf['specimen_date'])+1-minday
for (date,region,var,n) in zip(sgtf['specimen_date'],sgtf['UKHSA_region'],sgtf['sgtf'],sgtf['n']):
  day=Date(date)
  daynum=day-minday
  if daynum>=0 and daynum<nsgtf:
    place=region
    if place not in vocnum: vocnum[place]=np.zeros([nsgtf,2],dtype=int)
    vocnum[place][daynum][int("SGTF" in var)]+=n
  if date>=Date('2021-10-01') and date<Date('2021-11-10'):
    background[int("SGTF" in var)]+=n

# Adjust for non-Omicron SGTFs, based on the assumption that these are in a non-location-dependent proportion to the number of non-Omicron cases
f=background[1]/background[0]
for place in vocnum:
  for daynum in range(nsgtf):
    vocnum[place][daynum][1]=max(vocnum[place][daynum][1]-int(f*vocnum[place][daynum][0]+.5),0)

vocnum['England']=sum(vocnum.values())
nsgtf=vocnum['England'].shape[0]

step=7
with open('gg-by-region%d'%step,'w') as fp:
  n=maxday-minday
  #n=min(nsgtf,nspec)
  for loc in vocnum:
    if loc=='England': continue
    vv=vocnum[loc][:n,:]
    cv=casesbyregion[loc][:n,:].sum(axis=1)[:,None]*(vv/vv.sum(axis=1)[:,None])
    #cv=vv
    for i in range(n-step):
      if (vv[i:i+2*step:step,:]>=50).all():
        print("%7.4f %7.4f"%(log(cv[i+step][0]/cv[i][0])/step,log(cv[i+step][1]/cv[i][1])/step),Date(minday+i),Date(minday+i+step),"%6d %6d %6d %6d"%tuple(vv[i:i+2,:].reshape(-1)),loc,file=fp)
