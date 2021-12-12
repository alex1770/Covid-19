from stuff import *
from math import exp,log,sqrt
from scipy.optimize import minimize
import numpy as np
from random import random, seed

cases=loadcsv('SAcasecounts.csv')
cases['South Africa']=cases['Total']
del cases['Total']

mindate='2021-10-01'
minday=datetoday(mindate)
day0=datetoday(cases['Date'][0])
if minday>day0:
  for x in cases: cases[x]=cases[x][minday-day0:]
N=len(cases['Date'])
day0=datetoday(cases['Date'][0])
day1=datetoday('2021-11-01')
monday=datetoday('2021-11-01')
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
xx=[5, -0.026, 0, 0.22, 0,0,0,0,0,0,0]
bounds=[(xx[0]-10,xx[0]+10), (-0.5,0.1), (xx[2]-10,xx[2]+10), (0.05, 0.8), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (0,0)]

# Optimise after perturbing target data (cases)
def opt(pert=0):
  targ=[log(x+eps)+sqrt(pert/(x+eps))*(random()*2-1) for x in cases[locname]]
  res=minimize(NLL,xx,args=(targ,),bounds=bounds,method="SLSQP",options={'maxiter':10000})#, 'eps':1e-4, 'ftol':1e-12})
  if not res.success: raise RuntimeError(res.message)
  return res.x,res.fun

# Central estimate
params0,f0=opt(0)
l=expand(params0)
targ=[log(x+eps) for x in cases[locname]]
#for d in range(N):
#  print(daytodate(day0+d),"%8.1f  %8.1f  %6.3f  %6.3f   %6.3f"%(l[d][0],l[d][1],l[d][3],targ[d],l[d][3]-targ[d]))
#print()

# Work out residuals to see how 'accurate' the input data is.
# resid measures how 'pseudo' the pseudo-Poisson residuals are.
resid=sum((l[d][3]-targ[d])**2*exp(targ[d]) for d in range(N))/N
perturbation=resid*3# Scale up so that the perturbuation distribution U[-sqrt(perturbation),sqrt(perturbation)] has variance equal to resid

# Reset first guess to central estimate
xx=params0
bounds=[(xx[0]-5,xx[0]+5), (xx[1]-0.3,xx[1]+0.3), (xx[2]-5,xx[2]+5), (xx[3]-0.3,xx[3]+0.3), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (-1,1), (0,0)]

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

pr("Growth advantage of Omicron over Delta: %.0f%% (%.0f%% - %.0f%%) each day = continuous growth of %.2f (%.2f - %.2f) per day"%((exp(params0[3]-params0[1])-1)*100,(exp(growthadv[0])-1)*100,(exp(growthadv[1])-1)*100,params0[3]-params0[1],growthadv[0],growthadv[1]))
pr("Approximate R_t(Omicron)/R_t(Delta): %.2f (%.2f - %.2f), assuming a generation time of %.1f days"%(exp((params0[3]-params0[1])*gentime),exp(growthadv[0]*gentime),exp(growthadv[1]*gentime),gentime))
#print('"'+locname+'"',exp((params0[3]-params0[1])*gentime),exp(growthadv[0]*gentime),exp(growthadv[1]*gentime),file=sys.stderr)
#print('"'+locname+'"',params0[1],params0[3],file=sys.stderr)
pr("Doubling time of Omicron/Delta: %.1f (%.1f - %.1f) days"%(log(2)/(params0[3]-params0[1]),log(2)/growthadv[1],log(2)/growthadv[0]))
pr("")
pr("Growth of Omicron: %.0f%% (%.0f%% - %.0f%%) each day = continuous growth of %.2f (%.2f - %.2f) per day"%((exp(params0[3])-1)*100,(exp(growth[0])-1)*100,(exp(growth[1])-1)*100,params0[3],growth[0],growth[1]))
pr("Approximate R_t(Omicron): %.2f (%.2f - %.2f), assuming a generation time of %.1f days"%(exp(params0[3]*gentime),exp(growth[0]*gentime),exp(growth[1]*gentime),gentime))
pr("Doubling time of Omicron: %.1f (%.1f - %.1f) days"%(log(2)/params0[3],log(2)/growth[1],log(2)/growth[0]))
pr("")
pr("Description: http://sonorouschocolate.com/covid19/index.php?title=Early\\\_Omicron\\\_Growth\\\_Estimate")
pr("Data source: https://www.nicd.ac.za/ at "+cases['Date'][-1])

# Make graphs
os.makedirs(outputdir,exist_ok=True)
l=expand(params0,N)
for adj in [0,1]:
  if adj: weekadj=np.exp(params0[4:]-sum(params0[4:])/7)
  else: weekadj=[1]*7
  
  title='Covid-19 case count in '+locname+' (with weekday adjustment)'*adj+' decomposed as sum of falling exponential for Delta + rising exponential for Omicron'
  data=[]
  data.append({
    'title': 'Reported case count',
    'values': [(daytodate(day0+d),cases[locname][d]/weekadj[(d+day0-monday)%7]) for d in range(N)],
    'with': ('points',1),
    'extra': 'pt 5'
  })
  data.append({
    'title': 'Model',
    'values': [(daytodate(day0+d),exp(l[d][3])/weekadj[(d+day0-monday)%7]) for d in range(N)]
  })
  label='set label at graph 0.25,0.98 "'+'\\n'.join(text)+'"'
  outfn=os.path.join(outputdir,locname.replace(' ','_')+'_cases'+'_adj'*adj+'.png')
  makegraph(title=title, data=data, ylabel='New cases per day', outfn=outfn, extra=[label,'set key left','set logscale y 2'])
