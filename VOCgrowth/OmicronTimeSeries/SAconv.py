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
if len(sys.argv)>1: locname=sys.argv[1].replace('_',' ')
if len(sys.argv)>2: ntrials=int(sys.argv[2])

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
while data[key][day4-dataday0]=='': day4=day4-1

target=data[key][day3-dataday0:day4-dataday0+1]

est_delta=[]
est_omicron=[]
weekdayadj=[]
for day in Daterange(day3,day4+1):
  weekdayadj.append(weekdayadj7[(day-monday)%7])
  if day<day2:
    est_delta.append(data['cases'][day-dataday0])
    est_omicron.append(0)
  else:
    ncases=data['cases'][day-dataday0]
    modelleddelta=exp(xx[0]+xx[1]*(day-day1))*weekdayadj[-1]
    est_delta.append(min(ncases,modelleddelta))
    est_omicron.append(ncases-est_delta[-1])

n=day4-day3+1
est_delta=np.array(est_delta)
est_omicron=np.array(est_omicron)
weekdayadj=np.array(weekdayadj)

# Make graphs
os.makedirs(outputdir,exist_ok=True)
offset=day0-day3
for adj in [0,1]:
  if adj: weekadj=weekdayadj[offset:]
  else: weekadj=[1]*N
  
  title='Covid-19 case count in '+locname+' (with weekday adjustment)'*adj+' decomposed as sum of falling exponential for Delta + rising exponential for Omicron'
  grdata=[]
  grdata.append({
    'title': 'Reported case count',
    'values': [(str(day0+d),expdat[d]/weekadj[d]) for d in range(N)],
    'with': ('points',1),
    'extra': 'pt 5'
  })
  grdata.append({
    'title': 'Model',
    'values': [(str(day0+d),(est_delta[offset+d]+0*est_omicron[offset+d])/weekadj[d]) for d in range(N)]
  })
  grdata.append({
    'title': 'Target',
    'values': [(str(day0+d),target[offset+d]) for d in range(N)]
  })
  label='set label at graph 0.25,0.98 "Stuff"'
  outfn=os.path.join(outputdir,locname.replace(' ','_')+'_cases'+'_adj'*adj+'.png')
  makegraph(title=title, data=grdata, ylabel='New cases per day', outfn=outfn, extra=[label,'set key left','set logscale y 2'])
