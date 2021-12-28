from stuff import *
from math import log,exp
from scipy.optimize import minimize
import numpy as np
np.set_printoptions(precision=3,linewidth=250,suppress=True)

l=[x for x in os.listdir('.') if x[:19]=='sgtf_regionepicurve']
if l==[]: raise RuntimeError("No sgtf_regionepicurve csv file found in current directory")
sgtf=loadcsv(max(l))

minday0=Date('2021-09-20')
minday=Date('2021-10-17')
pubday=Date('2021-12-27')
maxday=max(Date(d) for d in sgtf['specimen_date'])+1
ndays=maxday-minday
skip=2

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

if 0:
  for d in range(ndays):
    ndelta,nomicron=vocnum['England'][d]
    date=minday+d
    if date>='2021-11-15':
      print(date,'%6d %6d   %8.5f'%(ndelta,nomicron,nomicron/(ndelta+nomicron)))
  poi

location='London'
#ages=[(a,a+10) for a in range(0,70,10)]+[(70,150)]
#ages=[(0,40),(40,60),(60,150)]
ages=[(0,50),(50,150)]
#ages=[(20,40),(60,150)]
nages=len(ages)
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
vn=vocnum[location]
nsgtf=vn.shape[0]

#for day in Daterange(minday,maxday): print(day,"%8.3f"%(log(vn[day-minday][1]/vn[day-minday][0]+1e-100)),cases[day-minday])
#for day in Daterange(minday,maxday): print(day,cases[day-minday].sum()/vn[day-minday].sum())
#for day in Daterange(minday,maxday): print(day,sp[day-minday0].sum()/vn[day-minday].sum())

# Variables:
# nspec*nages: log(ecases) array of shape (nspec,nages)
# nages:       crossover days relative to minday
# 1:           h

# Initial guess
lcases0=np.log(cases+0.5)
cross0=[Date('2021-12-15')-minday]*nages
h0=0.3
xx0=np.concatenate([lcases0.reshape([-1]),cross0,[h0]])

flcases=lcases0.reshape([-1])
change=3
bounds0=list(zip(flcases-change,flcases+change))
bounds1=[(Date('2021-12-01')-minday,Date('2022-01-01')-minday)]*nages
bounds2=[(0.1,0.5)]
bounds=bounds0+bounds1+bounds2

assert (np.array(bounds)[:,0]<=xx0).all() and (xx0<=np.array(bounds)[:,1]).all()

nif1=0.3# Non-independence factor for case counts
nif2=1# Non-independence factor for SGTF counts
sfmin=50 # Penalty for growth decreasing
sfmax=2000# Penalty for growth increasing

# Axes: (day, age, variant)   variant=0 or 1
#    E.g., 36 x 3 x 2

def unpack(xx):
  # log of expected number of cases, summed over variant
  lcases=xx[:nspec*nages].reshape(lcases0.shape)

  # Crossover days for each age
  cross=xx[nspec*nages:nspec*nages+nages]

  # Variant growth advantage
  h=xx[nspec*nages+nages]

  return lcases,cross,h

def LL(xx):
  ll=0
  
  lcases,cross,h=unpack(xx)
  
  ecases=np.exp(lcases)
  
  # Scaled Poisson for cases
  ll+=((cases-ecases).sum()+(cases*(lcases-lcases0)).sum())*nif1
  
  # May be better to normalise like this, but requires using a constraint to "gauge-fix"
  # ecases*=cases.sum()/ecases.sum()
  # ecases.sum() is a constant; subtract off lcases0 to keep numbers not too big
  # ll+=(cases*(lcases-lcases0)).sum()*nif1
  
  # Variant odds ratios for each (day, age band)
  vf=np.exp((np.arange(nspec)[:,None]-np.array(cross)[None,:])*h)
  
  # Expected no. of cases for each variant
  v0=ecases*1/(1+vf)
  v1=ecases*vf/(1+vf)
  # Sum over ages
  s0=v0.sum(axis=1)
  s1=v1.sum(axis=1)
  
  # Binomial probability of variant
  pp=(s1/(s0+s1))[:nsgtf]
  ll+=(vn[:,0]*np.log(1-pp)+vn[:,1]*np.log(pp)).sum()*nif2
  
  gr=np.log(v0[1:]/v0[:-1])
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
print(res.message)
print('Function value =',res.fun)
print()
xx=res.x
for i in range(len(xx)):
  if xx[i]>bounds[i][1]-1e-3 or xx[i]<bounds[i][0]+1e-3: print("Variable %d = %g hitting bound (%g, %g)"%(i,xx[i],bounds[i][0],bounds[i][1]))

lcases,cross,h=unpack(xx)
ecases=np.exp(lcases)
vf=np.exp((np.arange(nspec)[:,None]-np.array(cross)[None,:])*h)
v0=ecases*1/(1+vf)
v1=ecases*vf/(1+vf)
s0=v0.sum(axis=1)
s1=v1.sum(axis=1)
gr=np.log(v0[1:]/v0[:-1])

print('Location:',location)
print('Omicron/Delta growth = %.3f'%h)
print('Age bands:',ages)
print("Crossovers:",[Date(minday+int(x+.5)) for x in cross])
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
  print(Date(minday+d),end='')
  for a in range(nages): print('     %6.0f %6.0f'%(v0[d][a],v1[d][a]),end='')
  print('  ',end='')
  for a in range(nages):
    if d<nspec-1:
      print('  %6.3f'%gr[d][a],end='')
      print('  %6.3f'%(log(v0[d][a]*v1[d+1][a]/(v0[d+1][a]*v1[d][a]))),end='')
    else:
      print('       -',end='')
      print('       -',end='')
  if d<nspec-1:
    print('  %6.3f'%(log(s0[d]*s1[d+1]/(s0[d+1]*s1[d]))),end='')
  else:
    print('       -',end='')
  print('  %7.3f'%log(s1[d]/s0[d]),end='')
  if d<nsgtf and vn[d][0]>0 and vn[d][1]>0: print('  %7.3f'%log(vn[d][1]/vn[d][0]),end='')
  else: print('        -',end='')
  print()
