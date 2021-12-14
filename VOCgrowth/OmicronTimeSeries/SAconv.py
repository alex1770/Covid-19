from stuff import *
from math import exp,log,sqrt
from scipy.optimize import minimize
import numpy as np
from random import random, seed

data=loadcsv('Gauteng.Samizdat.csv')

twoexpdates=[Date('2021-09-01'), Date('2021-12-04')]
monday=Date('2021-11-01')
day0=twoexpdates[0]
day1=Date('2021-11-03')# Centre exponentials on approx minimum for stability

expdat=[]
for (date,numc) in zip(data['date'],data['cases']):
  if date>=twoexpdates[0] and date<=twoexpdates[1]: expdat.append(numc)
N=len(expdat)
#expdat=[100*exp(-0.05*(d-day1))+200*exp(0.25*(d-day1)) for d in Daterange(twoexpdates[0],twoexpdates[1]+1)]
  
np.set_printoptions(precision=4,linewidth=1000)
locname='South Africa'
outputdir='output'
conf=0.95
ntrials=1000
eps=0.1# Add this virtual case to case count to prevent singularity
seed(42)
if len(sys.argv)>1: locname=sys.argv[1]
if len(sys.argv)>2: ntrials=int(sys.argv[2])
locname=locname.replace('_',' ')
locname1=locname.replace(' ','_')

# Model
#         Delta                 Omicron
# log(exp(x0+x1*(day-day1))+exp(x2+x3*(day-day1))) + x_{4+(day-monday)%7}
# x1<0, x3>0
# day1 = 2021-11-01 (approximate minimum of cases)
# x0 ~= 4
# x2 ~= 2
# x10  = 0 (gauge fix)

def expand(xx,n=N):
  l=[]
  for d in range(n):
    n_delta=exp(xx[0]+xx[1]*(d+day0-day1)+xx[4+(day0+d-monday)%7])
    n_omicron=exp(xx[2]+xx[3]*(d+day0-day1)+xx[4+(day0+d-monday)%7])
    logn=log(n_delta+n_omicron)
    l.append((n_delta,n_omicron,xx[4+(day0+d-monday)%7],logn))
  return l

def NLL(xx,targ):
  l=expand(xx)
  s=0
  for d in range(N):
    # Don't try to match in the crossover period, because Omicron is going to be "stochastic" there
    s+=(l[d][3]-targ[d])**2/2*(abs(day0+d-day1)>=10)
    #s+=(l[d][3]-targ[d])**2/2*abs(day0+d-day1)
    #s+=(l[d][3]-targ[d])**2/2*(targ[d]/targ[day1-day0])
    #s+=(exp(l[d][3])-exp(targ[d])*(l[d][3]+1-targ[d]))/1000
    #s+=(l[d][3]-targ[d])**2/2
  return s

xx=[5, -0.026, 0, 0.22, 0,0,0,0,0,0,0]
bounds=[(xx[0]-10,xx[0]+10), (-0.5,0.1), (xx[2]-10,xx[2]+10), (0.05, 0.8), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (0,0)]

# Central estimate
targ=[log(x+eps) for x in expdat]
res=minimize(NLL,xx,args=(targ,),bounds=bounds,method="SLSQP",options={'maxiter':10000})#, 'eps':1e-4, 'ftol':1e-12})
if not res.success: raise RuntimeError(res.message)
xx,f0=res.x,res.fun
adj=log(sum(expdat[d]/exp(xx[4+(day0+d-monday)%7]) for d in range(N))/sum(expdat))
weekdayadj7=[exp(xx[4+r]+adj) for r in range(7)]

# Now going to wider date range to infer hospitalisation and death dependence on cases
# We're going to assume anything before day2 (2021-10-18) is non-Omicron
# After that, we're going to use extrapolate the exponential fall +weekday adjustment for non-Omicron (assumed Delta)
# and then whatever remains is assumed to be Omicron
day2=Date('2021-10-18')
day3=Date('2021-01-01')
day4=Date(max(data['date']))
dataday0=data['date'][0]
key='admissions'
while data[key][day4-dataday0]=='': day4-=1

target=[float(x) for x in data[key][day3-dataday0:day4-dataday0+1]]
# We're not going to predict some initial segment of target because of convolution width

delta=[]
omicron=[]
weekdayadj=[]
for day in Daterange(day3,day4+1):
  weekdayadj.append(weekdayadj7[(day-monday)%7])
  if day<day2:
    delta.append(data['cases'][day-dataday0])
    omicron.append(0)
  else:
    ncases=data['cases'][day-dataday0]
    modelleddelta=exp(xx[0]+xx[1]*(day-day1))*weekdayadj[-1]
    delta.append(min(ncases,modelleddelta))
    omicron.append(ncases-delta[-1])

n=day4-day3+1
delta=np.array(delta)
omicron=np.array(omicron)
weekdayadj=np.array(weekdayadj)

title='Estimated new cases per day of non-Omicron and Omicron in '+locname
minc=10
grdata=[]
grdata.append({
  'title': 'Estimated non-Omicron',
  'values': [(str(day3+d),delta[d]) for d in range(n) if delta[d]>=minc]
})
grdata.append({
  'title': 'Estimated Omicron',
  'values': [(str(day3+d),omicron[d]) for d in range(n) if omicron[d]>=minc]
})
#label='set label at graph 0.25,0.98 "Stuff"'
outfn=os.path.join(outputdir,locname1+'_estimatedvariants.png')
makegraph(title=title, data=grdata, ylabel='New cases per day', outfn=outfn, extra=['set key left', 'set logscale y 2'])


from scipy.stats import gamma
def getkernel(shape, scale, length):
  tt=[]
  for i in range(length+1):
    x=gamma.sf(i,shape,scale=scale)
    tt.append(x)
    #if x<cutoff: break
  tt=np.array(tt)
  kernel=tt[:-1]-tt[1:]
  return kernel/kernel.sum()

def getprobs(cdelta,comicron,target,alpha):
  m=len(cdelta)
  assert len(comicron)==m and len(target)==m
  # Coding of 0,1,...,2m-1 is 0,...,m-1 <-> pdelta[], and m,...,2m-1 <-> pomicron
  A=np.zeros([2*m,2*m])
  b=np.zeros(2*m)
  for i in range(m):
    A[i,i]+=cdelta[i]**2
    A[m+i,m+i]+=comicron[i]**2
    A[i,m+i]+=cdelta[i]*comicron[i]
    A[m+i,i]+=cdelta[i]*comicron[i]
    b[i]+=cdelta[i]*target[i]
    b[m+i]+=comicron[i]*target[i]
  for i in range(m-1):
    A[i,i]+=alpha
    A[i+1,i+1]+=alpha
    A[i,i+1]-=alpha
    A[i+1,i]-=alpha
    A[m+i,m+i]+=alpha
    A[m+i+1,m+i+1]+=alpha
    A[m+i,m+i+1]-=alpha
    A[m+i+1,m+i]-=alpha
  pp=np.linalg.lstsq(A,b)[0]
  pdelta=pp[:m]
  pomicron=pp[m:]
  resid=((pdelta*cdelta+pomicron*comicron-target)**2).sum()
  return pdelta,pomicron,resid

alpha=1e9# Controls timescale on which probs can change
r=50# Max size of kernel

for shape in np.arange(0.5,2,0.25):
  for scale in np.arange(3,8,0.25):
    k=getkernel(shape,scale,r+1)
    assert r==len(k)-1
    cdelta=np.convolve(delta,k,'valid')
    comicron=np.convolve(omicron,k,'valid')
    pdelta,pomicron,resid=getprobs(cdelta,comicron,target[r:],alpha)
    print("%6.3f %6.3f  %12g"%(shape,scale,resid))

title='Probs'
grdata=[]
grdata.append({
  'title': 'P(hosp|non-omicron case)',
  'values': [(str(day3+d+r),pdelta[d]) for d in range(n-r)]
})
grdata.append({
  'title': 'P(hosp|omicron case)',
  'values': [(str(day3+d+r),pomicron[d]) for d in range(n-r)]
})
#label='set label at graph 0.25,0.98 "Stuff"'
outfn=os.path.join(outputdir,locname1+'_probs.png')
makegraph(title=title, data=grdata, ylabel='P(hosp|case)', outfn=outfn, extra=['set key left'])

# Make graphs
os.makedirs(outputdir,exist_ok=True)
offset=day0-day3

if 1:
  title='Thing'
  grdata=[]
  grdata.append({
    'title': 'Delta',
    'values': [(str(day3+d),delta[d]) for d in range(n)]
  })
  grdata.append({
    'title': 'Omicron',
    'values': [(str(day3+d),omicron[d]) for d in range(n)]
  })
  grdata.append({
    'title': 'Target',
    'values': [(str(day3+r+d),target[r+d]) for d in range(n-r)],
    #'with': ('points',1),
    #'extra': 'pt 5'
  })
  grdata.append({
    'title': 'Prediction',
    'values': [(str(day3+r+d),cdelta[d]) for d in range(n-r)]
  })
  label='set label at graph 0.25,0.98 "Stuff"'
  outfn=os.path.join(outputdir,locname1+'_predictor.png')
  makegraph(title=title, data=grdata, ylabel='New cases per day', outfn=outfn, extra=[label,'set key left','set logscale y 2'])
