from stuff import *
from math import log,exp
from scipy.optimize import minimize
from random import randrange
from math import sqrt
import numpy as np
np.set_printoptions(precision=6,linewidth=250,suppress=True)

l=[x for x in os.listdir('.') if x[:19]=='sgtf_regionepicurve']
if l==[]: raise RuntimeError("No sgtf_regionepicurve*.csv file found in current directory; download from https://www.gov.uk/government/publications/covid-19-omicron-daily-overview")
sgtf=loadcsv(max(l))
regions=sorted(list(set(sgtf['UKHSA_region'])))
minday=Date('2021-11-25')
maxday=Date('2021-12-25')# Only go up to dates strictly before this one
pubday=getpublishdate()
discard=3# Discard last few case counts by specimen date since these are incomplete (irrelevant here because we're stopping much earlier anyway)
mincount=50
step=7

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

data=[]
wt=[]
with open('gg-by-region%d'%step,'w') as fp:
  n=maxday-minday
  #n=min(nsgtf,nspec)
  for loc in vocnum:
    if loc=='England': continue
    vv=vocnum[loc][:n,:]
    cv=casesbyregion[loc][:n,:].sum(axis=1)[:,None]*(vv/vv.sum(axis=1)[:,None])
    #cv=vv
    for i in range(n-step):
      if (vv[i:i+2*step:step,:]>=mincount).all():
        gr0=log(cv[i+step][0]/cv[i][0])/step
        gr1=log(cv[i+step][1]/cv[i][1])/step
        wt0=1/(1/vv[i,0]+1/vv[i+step,0])
        wt1=1/(1/vv[i,1]+1/vv[i+step,1])
        data.append((gr0,gr1))
        wt.append((wt0,wt1))
        print("%7.4f %7.4f"%(data[-1]),Date(minday+i),Date(minday+i+step),"%6d %6d %6d %6d"%tuple(vv[i:i+2*step:step,:].reshape(-1)),loc,file=fp)
data=np.array(data)
wt=np.array(wt)

def regress(data,wt):
  # Start with simple weighted regression, pretending uncertainty in (x,y) is isotropic, to give a starting point
  W=wt.sum(axis=1)
  (X,Y)=data[:,0],data[:,1]
  m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
  r=np.array([sum(W*Y),sum(W*X*Y)])
  c=np.linalg.solve(m,r)
  #print(c)
  #res=c[0]+c[1]*X-Y
  #tres=(W*res*res).sum()

  # Calculate weighted square-distance to putative best fit line of gradient b
  def f(b):
    (X,Y)=data[:,0],data[:,1]
    (V,W)=wt[:,0],wt[:,1]
    N=Y-b*X
    D=W+b*b*V
    a=sum(N/D)/sum(1/D)
    return (a,sum((N-a)*(N-a)/(2*D)))
  
  f0=f(c[1])[1]
  def g(b): return f(b)[1]/f0
  
  bounds=(c[1]*0.5,c[1]/0.5)
  res=minimize(g,[c[1]],bounds=[bounds],method="SLSQP")
  if not res.success: raise RuntimeError(res.message)
  b=res.x[0]
  (a,s)=f(b)
  #if b<bounds[0]+1e-6 or b>bounds[1]-1e-6: raise RuntimeError("Hit bounds %g %g"%bounds)
  #if b<bounds[0]+1e-6 or b>bounds[1]-1e-6: print("Warning: Hit bounds %g %g"%bounds,file=sys.stderr)
  return (a,b)
  
# Bootstrap to get confidence intervals
n=len(data)
nsamp=10000
conf=0.95

samples=[]
print("Generating %d samples"%nsamp)
for samp in range(nsamp):
  data1=np.zeros(data.shape)
  wt1=np.zeros(wt.shape)
  for i in range(n):
    r=randrange(n)
    data1[i]=data[r]
    wt1[i]=wt[r]
  samples.append(regress(data1,wt1))

# Simple CI of gradient
samples.sort(key=lambda x:x[1])
low=samples[int((1-conf)/2*nsamp)]
high=samples[int((1+conf)/2*nsamp)]
central=regress(data,wt)
#print(central)
#print(low,high)
#print()

# Compare with Delta

# https://www.medrxiv.org/content/10.1101/2021.10.21.21265216v1
# Delta (intrinsic): 4.6 (4.0-5.4) days  3.1 (3.0-3.7) days
# >>> from scipy.stats import gamma
# >>> a=170;mu=4.64;gamma.ppf(.025,a,scale=mu/a),gamma.ppf(.975,a,scale=mu/a)
# (3.9687044779630987, 5.3629770162251926)

from scipy.stats import gamma,expon
csamples=[]
for (a,b) in samples:
  # Generate meanGT_Delta with characteristics 4.6 (4.0 - 5.4)
  ag=170;mu=4.64;mean_d=gamma.rvs(ag,scale=mu/ag)

  # Generate sdGT_Delta with characteristics 3.1 (3.0 - 3.7)
  # This is very skew, so do it in two pieces, setting the median (not the mean) to 3.1
  e=expon.ppf(0.95)
  if randrange(2): sd_d=3.1+.6/e*expon.rvs()
  else: sd_d=3.1-.1/e*expon.rvs()
  
  th_d=sd_d**2/mean_d
  k_d=(mean_d/sd_d)**2
  th_o=1/(b/th_d-a)
  k_o=k_d# Assume equal shape (or coefficient of variation)
  rr=(b*th_o/th_d)**k_d
  mean_o=th_o*k_o
  sd_o=th_o*sqrt(k_o)
  gtr=mean_o/mean_d
  csamples.append([mean_d,sd_d,mean_o,sd_o,gtr,rr])
csamples=np.array(csamples)

desc=["meanGT_Delta","sdGT_Delta","meanGT_Omicron","sdGT_Omicron","GT_Om/GT_Delta","R_Om/R_Delta"]
nd=len(desc)
for i in range(nd):
  lsamp=list(csamples[:,i])
  lsamp.sort()
  low=lsamp[int((1-conf)/2*nsamp)]
  med=lsamp[int(0.5*nsamp)]
  high=lsamp[int((1+conf)/2*nsamp)]
  print("%-16s  %6.3f (%.3f - %.3f)"%(desc[i],med,low,high))
  
