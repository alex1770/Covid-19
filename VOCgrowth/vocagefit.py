from stuff import *
from math import log,exp
from scipy.optimize import minimize
import numpy as np
np.set_printoptions(precision=3,linewidth=250,suppress=True)

l=[x for x in os.listdir('.') if x[:19]=='sgtf_regionepicurve']
if l==[]: raise RuntimeError("No sgtf_regionepicurve*.csv file found in current directory; download from https://www.gov.uk/government/publications/covid-19-omicron-daily-overview")
sgtf=loadcsv(max(l))

location='England'
#ages=[(a,a+10) for a in range(0,70,10)]+[(70,150)]
#ages=[(0,40),(40,60),(60,150)]
ages=[(0,50),(50,150)]
#ages=[(0,60),(60,70),(70,150)]
#ages=[(20,40),(60,150)]
nages=len(ages)

minday0=Date('2021-09-20')
minday=Date('2021-10-17')
pubday=getpublishdate()
maxday=max(Date(d) for d in sgtf['specimen_date'])+1
ndays=maxday-minday
skip=2
useagedata=True

sp0,sp=getcasesbyagespeccomplete(minday=minday0,maxday=pubday,ages=ages,location=location)

# Simple weekday adjustment by dividing by the average count for that day of the week.
# Use a relatively stable period (inclusive) over which to take the weekday averages.
weekadjdates=['2021-09-20','2021-11-21']
weekadj=np.zeros(7)
for day in Daterange(*weekadjdates):
  assert day>=minday0 and day<pubday-7
  weekadj[day%7]+=sp0[day-minday0].sum()
weekadjp=weekadj*7/sum(weekadj)
cases=np.array([sp[day-minday0]/weekadjp[day%7] for day in Daterange(minday,pubday-skip)])
nspec=cases.shape[0]

# Get SGTF data into a suitable form
vocnum={}
background=[0,0]
for (date,region,var,n) in zip(sgtf['specimen_date'],sgtf['UKHSA_region'],sgtf['sgtf'],sgtf['n']):
  day=Date(date)
  daynum=day-minday
  if daynum>=0 and daynum<ndays:
    place=region
    if place not in vocnum: vocnum[place]=np.zeros([ndays,2],dtype=int)
    vocnum[place][daynum][int("SGTF" in var)]+=n
  if date>=Date('2021-10-01') and date<Date('2021-11-10'):
    background[int("SGTF" in var)]+=n
# Adjust for non-Omicron SGTFs, based on the assumption that these are in a non-location-dependent proportion to the number of non-Omicron cases
f=background[1]/background[0]
for place in vocnum:
  for daynum in range(ndays):
    vocnum[place][daynum][1]=max(vocnum[place][daynum][1]-int(f*vocnum[place][daynum][0]+.5),0)
vocnum['England']=sum(vocnum.values())
vn=vocnum[location]
nsgtf=vn.shape[0]

# From fig. 8B of Tech Briefing 33, https://www.gov.uk/government/publications/investigation-of-sars-cov-2-variants-technical-briefings
avd=loadcsv('age_variant_date.csv')
navd=max(Date(x+'US') for x in avd['spec.date'])+1-minday
avdata=np.zeros([navd,nages,2],dtype=int)
for (date,var,age,n) in zip(avd['spec.date'],avd['variant'],avd['age.group'],avd['count_age']):
  a=parseage(age)
  l=[i for (i,b) in enumerate(ages) if a[0]>=b[0] and a[1]<=b[1]]
  assert len(l)==1
  d=Date(date+'US')-minday
  if d>=0:
    avdata[d,l[0],int(var=='Omicron')]+=n
if useagedata: print("Using age-variant data")
else: print("Not using age-variant data")

# Normalise variant-age data to sgtf data with which it disagrees because I think the sgtf data is more reliable
# But keep overall normalisation so that the optimisation process correctly weights according to how much evidence there is
# I.e., (for all dates) make sum_{age}avdata(age,variant) = const*sgtfdata(variant)
# where const is chosen so that sum_{age,variant}avdata(age,variant) is unchanged.
a=avdata.astype(float)+1e-20
s=vn[:navd]+1e-20
avdata=a/a.sum(axis=1)[:,None,:]*a.sum(axis=(1,2))[:,None,None]*s[:,None,:]/s.sum(axis=1)[:,None,None]

if 0:
  for d in range(ndays):
    ndelta,nomicron=vocnum['England'][d]
    date=minday+d
    if date>='2021-11-15':
      print(date,'%6d %6d   %8.5f'%(ndelta,nomicron,nomicron/(ndelta+nomicron)))
  poi


with open('gg-by-region','w') as fp:
  n=min(nsgtf,nspec)
  for loc in vocnum:
    if loc=='England': continue
    vv=vocnum[loc][:n,:]
    #cv=casesbyregion[:n,None]*(vv/vv.sum(axis=1)[:,None])
    cv=vv
    for i in range(n-1):
      if (vv[i:i+2,:]>=100).all():
        print("%7.4f %7.4f"%(log(cv[i+1][0]/cv[i][0]),log(cv[i+1][1]/cv[i][1])),file=fp)

if 0:
  with open('gg-by-age','w') as fp:
    assert location=='England'
    n=min(navd,nspec)
    vv=avdata[:n,:,:]# Date, age, variant
    cv=vv
    #cv=cases[:n,:,None]*(vv/vv.sum(axis=2)[:,:,None])
    for a in range(nages):
      for i in range(n-1):
        if (vv[i:i+2,a,:]>=10).all():
          print("%7.4f %7.4f"%(log(cv[i+1,a,0]/cv[i,a,0]),log(cv[i+1,a,1]/cv[i,a,1])),file=fp)

pivotdate=Date('2021-12-16')
p=pivotdate-minday

# Variables:
# nspec*nages: log(ecases_delta) (array of shape (nspec,nages))
# nages:       log(ecases_omicron) on pivotdate
# 1:           h

# Axes: (day, age, variant)   variant=0 or 1
#    E.g., 36 x 3 x 2

def unpack(xx):
  # lcases = log of expected number of cases
  # h = variant growth advantage
  lcases=np.zeros([nspec,nages,2])
  lcases[:,:,0]=xx[:nspec*nages].reshape([nspec,nages])
  lcasesop0=xx[nspec*nages:nspec*nages+nages]
  h=xx[nspec*nages+nages]
  rho=xx[nspec*nages+nages+1]
  lcases[:,:,1]=lcasesop0[None,:]+(lcases[:,:,0]-lcases[p,:,0]+np.arange(-p,nspec-p)[:,None]*h)/rho
  ecases=np.exp(lcases)
  return lcases,ecases,h,rho

# Initial guess (could improve this)
lcasesd0=np.log(cases+0.5);lcasesd0[p:,:]=lcasesd0[p,None]
lcasesop0=np.copy(lcasesd0[p,:])
h0=0.3
rho0=1
xx0=np.concatenate([lcasesd0.reshape([-1]),lcasesop0,[h0,rho0]])
lcases0,ecases0,h0,rho0=unpack(xx0)
esv0=ecases0.sum(axis=2)
lesv0=np.log(esv0)

change=10
bounds0=list(zip(xx0[:-2]-change,xx0[:-2]+change))
bounds1=[(0.1,0.5)]
bounds2=[(0.1,2)]
bounds=bounds0+bounds1+bounds2

assert (np.array(bounds)[:,0]<=xx0).all() and (xx0<=np.array(bounds)[:,1]).all()

nif1=0.3# Non-independence factor for case counts
nif2=1# Non-independence factor for SGTF counts
sfmin=100 # Penalty for growth decreasing
sfmax=2000 # Penalty for growth increasing

def LL(xx):
  ll=0
  
  lcases,ecases,h,rho=unpack(xx)
  
  # Scaled Poisson for cases
  esv=ecases.sum(axis=2)
  lesv=np.log(esv)
  ll+=((cases-esv).sum()+(cases*(lesv-lesv0)).sum())*nif1
  
  # Sum over ages
  aa=ecases.sum(axis=1)
  
  # Binomial probability of variant
  p1=(aa[:,1]/aa.sum(axis=1))[:nsgtf]
  ll+=(vn[:,0]*np.log(1-p1)+vn[:,1]*np.log(p1)).sum()*nif2

  if useagedata:
    # Binomial probability of variant, by age
    lp=np.log(ecases/ecases.sum(axis=2)[:,:,None])
    n=min(nspec,navd)
    ll+=(avdata[:n]*lp[:n]).sum()*nif2

  # Distribution on change in growth rates
  gr=lcases[1:,:,0]-lcases[:-1,:,0]
  dgr=gr[1:,:]-gr[:-1,:]
  dgrmax=np.maximum(dgr,0)
  dgrmin=np.minimum(dgr,0)
  ll+=-(sfmin**2/2*(dgrmin*dgrmin).sum()+sfmax**2/2*(dgrmax*dgrmax).sum())
  
  return ll

LL0=LL(xx0)
def NLL(xx):
  return (LL0-LL(xx))/1000

optmethod="SLSQP";minopts={"maxiter":10000}#,"eps":1e-4,'ftol':1e-12}
#optmethod="L-BFGS-B";minopts={"maxiter":10000,"maxfun":1000000}

res=minimize(NLL,xx0,bounds=bounds,method=optmethod,options=minopts)
if not res.success: raise RuntimeError(res.message)

if 0:
  # Redo
  xx0=res.x
  lcases0,ecases0,h0,rho0=unpack(xx0)
  esv0=ecases0.sum(axis=2)
  lesv0=np.log(esv0)
  LL0=LL(xx0)
  res=minimize(NLL,xx0,bounds=bounds,method=optmethod,options=minopts)
  if not res.success: raise RuntimeError(res.message)

print(res.message)
print('Function value =',res.fun)
print()
xx=res.x
for i in range(len(xx)):
  if xx[i]>bounds[i][1]-1e-3 or xx[i]<bounds[i][0]+1e-3: print("Variable %d = %g hitting bound (%g, %g)"%(i,xx[i],bounds[i][0],bounds[i][1]))

lcases,ecases,h,rho=unpack(xx)
aa=ecases.sum(axis=1)
gr=lcases[1:,:,0]-lcases[:-1,:,0]

def findcross(l):
  n=l.shape[0]
  for i in range(n):
    if l[i,1]>l[i,0]:
      if i==0: return 0
      return i-(l[i-1,1]-l[i-1,0]+l[i,1]-l[i,0]>0)
  return n-1

print('Location:',location)
print('Omicron/Delta growth = %.3f'%h)
print('Generation time ratio GT(Omicron)/GT(Delta) = %.3f'%rho)
print('Age bands:',ages)
print("Crossovers:",[minday+findcross(lcases[:,a,:]) for a in range(nages)])
print()

print("Date      ",end='')
for a in range(nages): print('       D_A%d   O_A%d'%(a,a),end='')
print('  ',end='')
for a in range(nages):
  print('   gr_A%d'%a,end='')
  print('   cr_A%d'%a,end='')
print('   cr_A*',end='')
print('  l(O/D)m',end='')
print('  l(O/D)a',end='')
print()

for d in range(nspec):
  print(minday+d,end='')
  for a in range(nages): print('     %6.0f %6.0f'%(ecases[d,a,0],ecases[d,a,1]),end='')
  print('  ',end='')
  for a in range(nages):
    if d<nspec-1:
      print('  %6.3f'%gr[d][a],end='')
      print('  %6.3f'%(lcases[d,a,0]+lcases[d+1,a,1]-lcases[d+1,a,0]-lcases[d,a,1]),end='')
    else:
      print('       -',end='')
      print('       -',end='')
  if d<nspec-1:
    print('  %6.3f'%(log(aa[d,0]*aa[d+1,1]/(aa[d+1,0]*aa[d,1]))),end='')
  else:
    print('       -',end='')
  print('  %7.3f'%log(aa[d,1]/aa[d,0]),end='')
  if d<nsgtf and vn[d][0]>0 and vn[d][1]>0: print('  %7.3f'%log(vn[d][1]/vn[d][0]),end='')
  else: print('        -',end='')
  print()

with open('gg-%s-modelled'%location,'w') as fp:
  n=min(nsgtf,nspec)
  for d in range(n-1):
    for a in range(nages):
      print('%7.4f %7.4f'%(lcases[d+1,a,0]-lcases[d,a,0],lcases[d+1,a,1]-lcases[d,a,1]),file=fp)

if location=='England':
  with open('av-comparison','w') as fp:
    n=min(nsgtf,nspec,navd)
    vv=avdata[:n,:,:]# Date, age, variant
    #cv=vv
    cv=cases[:n,:,None]*(vv/vv.sum(axis=2)[:,:,None])
    for d in range(n-1):
      if (vv[d:d+2,:,:]>=10).all():
        print(minday+d,end='',file=fp)
        for a in range(nages):
          for var in [0,1]:
            print('   %6.0f %6.0f'%(ecases[d,a,var],cv[d,a,var]),end='',file=fp)
        print(file=fp)
