from stuff import *
import json,requests,csv,os
import numpy as np
from math import log,exp
from scipy.optimize import minimize

# Ref https://twitter.com/PaulMainwood/status/1458860291492569101
#     https://www.gov.uk/government/publications/covid-19-vaccine-weekly-surveillance-reports (UKHSA)
#     Makes use of Paul's collation of UKHSA's data from https://github.com/PaulMainwood/vaccine-surveillance-reports

np.set_printoptions(precision=4,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=10000)

# ONS 2020 pop ages originally given as 0-18, 18-25, 25-30, ..., 75-80, 80+. Interpolate children age bands
# Get NIMSpop from dashboard (vaccinationsAgeDemographics metric)
ONSages=[(0,12),         (12,16),       (16,18),       (18,25), (25,30), (30,35), (35,40), (40,45), (45,50), (50,55), (55,60), (60,65), (65,70), (70,75), (75,80), (80,150)]
ONSpop= [12093288*12/18, 12093288*4/18, 12093288*2/18, 4709589, 3771493, 3824652, 3738209, 3476303, 3638639, 3875351, 3761782, 3196813, 2784300, 2814128, 2009992, 2855599]

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(10):
    try:
      response = requests.get(url+req, timeout=5)
      if response.ok: break
      error=response.text
    except BaseException as err:
      error=str(err)
  else: raise RuntimeError('Request failed: '+error)
  return response.json()['body'][::-1]

# Convert (eg) string ages '15-19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  x=x.strip()
  if x[-1]=='+': return (int(x[:-1]),150)
  if x[:6]=='Under ': return (0,int(x[6:]))
  x=x.replace('_to_','_').replace('-','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

# Output in muggle format
def unparseage(ages):
  (x,y)=ages
  if y<150: return "%d-%d"%(x,y-1)
  return "%d+"%x

#rawvax=get_data('areaType=nation&areaName=England&metric=vaccinationsAgeDemographics')

pmd=loadcsv("vaccine-surveillance-reports/vaccine surveillance.csv")
ukhsaages=sorted(parseage(a) for a in set(pmd['Age']))
graphdir='tempd'
os.makedirs(graphdir,exist_ok=True)

for ag in [(40,50)]:#ukhsaages:

  ty='Cases'
  #ty='Death 28 days'
  D0=[];D2=[]
  for (a,d,t,ul,d0,dhalf,d1,d2) in zip(pmd['Age'],pmd['Date'],pmd['Type'],pmd['Unlinked'],
                                       pmd['Not vaccinated'],
                                       pmd['Received one dose 1-20 days before specimen date'],
                                       pmd['Received one dose, 21 days or more before specimen date'],
                                       pmd['Second dose 14 days or more before specimen date']):
    if t!=ty: continue
    if parseage(a)!=ag: continue
    #print(a,d,t,ul,d0,dhalf,d1,d2)
    print(a,d,t,d0,d2)#,log(1-d2/d0)/log(2))
    D0.append(d0);D2.append(d2)
  
  # k=(vax pop)/(unvax pop)
  # Look for log(1-(1/k)d2(t)/d0(t)) to be linear
  # 1-(1/k)d2/d0 ~ exp(-alpha.(t-t0))*VE(t0)
  # d2/d0 ~ k*(1-exp(-alpha.(t-t0))*VE(t0))
  # d0 ~ B(d0+d2, 1/(1+k*(1-exp(-alpha.(t-t0))*VE(t0))))
  D0=np.array(D0);D2=np.array(D2)
  if D0.sum()==0 or D2.sum()==0: print();continue
  kmin=1.1*max(D2/D0)
  #kmin=max(kmin,10)#alter
  kmax=1000
  
  def NLL(xx):
    logk,al,VE0=xx
    LL=0
    k=exp(logk)
    for (t,(d0,d2)) in enumerate(zip(D0,D2)):
      p=1/(1+k*(1-exp(-al*t)*VE0))
      #print(k,al,VE0,t,p)
      p0=d0/(d0+d2)
      LL+=d0*log(p/p0)+d2*log((1-p)/(1-p0))
    return -LL/(sum(D0)+sum(D2))

  if 1:
    bounds=[(log(kmin),log(kmax)), (0,1), (0,0.9999)]
    res=minimize(NLL,(log(kmin*kmax)/2,0.06,0.7),bounds=bounds,method="SLSQP",options={'maxiter':10000})
    if not res.success: raise RuntimeError(res.message)
    print(res.message)
    logk,al,VE0=res.x
    k=exp(logk)
    print("k =",k)
    print("    P(vax) =",k/(1+k))
    print("Alpha =",al)
    print("VE(0) =",VE0)

  numgr=12
  ar=unparseage(ag)
  data=[]
  for i in range(numgr):
    k=exp(log(kmin)+(log(kmax)-log(kmin))*i/(numgr-1))
    values=[]
    for (t,(d0,d2)) in enumerate(zip(D0,D2)):
      dt=daytodate(datetoday('2021-08-19')+7*t)
      values.append((dt,k*log(1-d2/d0/k)))
    data.append({
      'title': 'k=%.3f'%k,
      'values': values,
    })
  title='Waning curve for ages %s'%ar
  makegraph(title=title, data=data, ylabel='k*log(1-d2/d0/k)', outfn=os.path.join(graphdir,'%s.png'%ar))
  #extra=["set key top left",'set style fill transparent solid 0.25'],ranges='[:] [0:1.1]', interval=86400*14)
  print()
