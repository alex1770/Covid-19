from stuff import *
from math import exp,log,sqrt
from scipy.optimize import minimize
import numpy as np
from random import random, seed

cases=loadcsv('SAcasecounts.csv')

N=len(cases['Date'])
day0=datetoday(cases['Date'][0])
day1=datetoday('2021-11-01')
monday=datetoday('2021-11-01')
np.set_printoptions(precision=4,linewidth=1000)
prov='Gauteng'
conf=0.95
ntrials=1000
seed(42)
if len(sys.argv)>1: prov=sys.argv[1]
if len(sys.argv)>2: ntrials=int(sys.argv[2])
locname=(prov if prov!='Total' else 'South Africa')

# Model
#         Delta                 Omicron
# log(exp(x0+x1*(day-day1))+exp(x2+x3*(day-day1))) + x_{4+(day-monday)%7}
# x1<0, x3>0
# day1 = 2021-11-01 (approximate minimum of cases)
# x0 ~= 5
# x2 ~= 2
# x10  = 0 (gauge fix)

def expand(xx,n=N):
  l=[]
  for d in range(n):
    n_delta=exp(xx[0]+xx[1]*(d+day0-day1))
    n_omicron=exp(xx[2]+xx[3]*(d+day0-day1))
    logn=log(n_delta+n_omicron)+xx[4+(day0+d-monday)%7]
    l.append((n_delta,n_omicron,xx[4+(day0+d-monday)%7],logn))
  return l

def NLL(xx,targ):
  l=expand(xx)
  s=0
  for d in range(N):
    s+=(l[d][3]-targ[d])**2/2
  return s

# First guess
xx=[5, -0.026, 2, 0.22, 0,0,0,0,0,0,0]
bounds=[(xx[0]-5,xx[0]+5), (-0.5,0.1), (xx[2]-5,xx[2]+5), (0.05, 0.4), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (0,0)]

# Optimise after perturbing target data (cases)
def opt(pert=0):
  targ=[log(x)+pert*(random()*2-1) for x in cases[prov]]
  res=minimize(NLL,xx,args=(targ,),bounds=bounds,method="SLSQP",options={'maxiter':10000})#, 'eps':1e-4, 'ftol':1e-12})
  if not res.success: raise RuntimeError(res.message)
  return res.x,res.fun

# Central estimate
params0,f0=opt(0)
l=expand(params0)
targ=[log(x) for x in cases[prov]]
#for d in range(N):
#  print(daytodate(day0+d),"%8.1f  %8.1f  %6.3f  %6.3f   %6.3f"%(l[d][0],l[d][1],l[d][3],targ[d],l[d][3]-targ[d]))
#print()

# Work out residuals to see how 'accurate' the input data is
resid=sqrt(sum((l[d][3]-targ[d])**2 for d in range(N))/N)
perturbation=resid*sqrt(3)# Scale up so that the perturbuation distribution U[-perturbation,perturbation] has standard deviation equal to resid

# Reset first guess to central estimate
xx=params0
bounds=[(xx[0]-5,xx[0]+5), (-0.5,0.1), (xx[2]-5,xx[2]+5), (0.05, 0.4), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (0,0)]

# Confidence analysis
params=[]
for i in range(ntrials):
  x0,f0=opt(perturbation)
  params.append(x0)
  #values.append(expand(x0))
  #print(log(2)/x0[3])

def getconf(l):
  l=sorted(l)
  n=len(l)
  def interp(x):
    i=int(x+.5)
    if i==0: return l[0]
    if i==n: return l[n-1]
    return (.5+i-x)*l[i-1]+(.5+x-i)*l[i]
  return interp((1-conf)/2*n),interp((1+conf)/2*n)

params=np.array(params)
growthadv=getconf(params[:,3]-params[:,1])
growth=getconf(params[:,3])
gentime=5
text=[]
def pr(x):
  print(x)
  text.append(x)

#pr("Confidence intervals generated using a perturbation of exp(-%.3f ... %.3f) to case counts"%(perturbation,perturbation))
pr("Growth advantage of Omicron over Delta: %.0f%% (%.0f%% - %.0f%%) per day"%((exp(params0[3]-params0[1])-1)*100,(exp(growthadv[0])-1)*100,(exp(growthadv[1])-1)*100))
pr("Approximate R_t(Omicron)/R_t(Delta): %.2f (%.2f - %.2f), assuming a generation time of %.1f days"%(exp((params0[3]-params0[1])*gentime),exp(growthadv[0]*gentime),exp(growthadv[1]*gentime),gentime))
pr("Growth of Omicron: %.0f%% (%.0f%% - %.0f%%) per day"%((exp(params0[3])-1)*100,(exp(growth[0])-1)*100,(exp(growth[1])-1)*100))
pr("Approximate R_t(Omicron): %.2f (%.2f - %.2f), assuming a generation time of %.1f days"%(exp(params0[3]*gentime),exp(growth[0]*gentime),exp(growth[1]*gentime),gentime))
pr("Doubling time of Omicron: %.1f (%.1f - %.1f) days"%(log(2)/params0[3],log(2)/growth[1],log(2)/growth[0]))
pr("Program used: https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/EarlyOmicronEstimate/twoexpestimate.py")
pr("Data source: https://www.nicd.ac.za/")

title='Case count in '+locname+' decomposed as sum of falling exponential for Delta + rising exponential for Omicron'
data=[]
data.append({
  'title': 'Actual case count',
  'values': [(daytodate(day0+d),log(cases[prov][d])/log(2)) for d in range(N)],
  'with': ('points',1),
  'extra': 'pt 5'
})
l=expand(params0,N)
data.append({
  'title': 'Model',
  'values': [(daytodate(day0+d),l[d][3]/log(2)) for d in range(N)]
})
label='set label at graph 0.25,0.98 "'+'\\n'.join(text)+'"'
makegraph(title=title, data=data, ylabel='log_2 new cases per day', outfn='SAcases.png', extra=[label,'set key left'])
