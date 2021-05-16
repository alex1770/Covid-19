from stuff import *
import sys
from scipy.optimize import minimize
from random import randrange,seed

# Get ltla_2021-05-15.csv from https://coronavirus.data.gov.uk/api/v2/data?areaType=ltla&metric=newCasesBySpecimenDate&format=csv

# Get SGTF/S-gene from fig. 13 from this
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/986378/Variants_of_Concern_Technical_Briefing_11_Data_England__1_.xlsx
# whose containing page is
# https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201

apicases=loadcsv("ltla_2021-05-15.csv")
sgtf=loadcsv("Variants_of_Concern_Technical_Briefing_11_Data_England-1-fig13.csv")

# Get case counts c0, c1, c2, c3 for weeks commencing on these dates
# Then estimate old R by (c1/c0)^(mgt/gap) and new R by (c3/c2)^(mgt/gap)
#weeks=['2021-04-11','2021-04-18', '2021-04-25','2021-05-02']
weeks=['2021-04-11','2021-04-18', '2021-04-28','2021-05-05']
#weeks=['2021-04-18','2021-04-25', '2021-04-28','2021-05-05']
mgt=5# Mean generation time in days
exclude=set()# set(["E08000001"])

days=[datetoday(w) for w in weeks]
cases={}
for (ltla,date,n) in zip(apicases['areaCode'],apicases['date'],apicases['newCasesBySpecimenDate']):
  if ltla not in cases: cases[ltla]=[0,0,0,0]
  day=datetoday(date)
  for i in range(4):
    if day>=days[i] and day<days[i]+7: cases[ltla][i]+=n

ltlas=set()
sgtfnum={}
for (ltla,n,r) in zip(sgtf['ltla.code'],sgtf['Number classifiable cases'],sgtf['Number S gene cases']):
  if n>0 and ltla in cases and ltla not in exclude:
    ltlas.add(ltla)
    sgtfnum[ltla]=(r,n)
ltlas=sorted(list(ltlas))

l=[]
with open('graph','w') as fp:
  for ltla in ltlas:
    (r,n)=sgtfnum[ltla]
    cc=cases[ltla]
    if cc[0]==0 or cc[1]==0 or cc[2]==0 or cc[3]==0: continue
    gapold=days[1]-days[0]
    gapnew=days[3]-days[2]
    Rold=(cc[1]/cc[0])**(mgt/gapold)
    Rnew=(cc[3]/cc[2])**(mgt/gapnew)
    x=r/n
    y=Rnew/Rold
    print("%5.3f %5.3f"%(x,y),file=fp)
    # Calculate approx variances of x and y
    p=(r+1)/(n+2);vx=p*(1-p)/n
    vy=y*y*sum(1/c for c in cc)*mgt*((1/gapold+1/gapnew)/2)# This mgt factor on the variance is a fudge
    l.append([x,y,vx,vy,n])

def err(xx,l,cutoff):
  a,b=xx
  tot=0
  for (x,y,vx,vy,n) in l:
    if n>=cutoff:
      tot+=(y-(a+b*x))**2/(vy+b*b*vx)
  return tot

for cutoff in [1,5,10,20,30,50]:
  xx=[1,1]
  bounds=[(0.5,1.5),(-1,2)]
  res=minimize(err,xx,args=(l,cutoff),method="SLSQP",bounds=bounds,options={"maxiter":1000})
  if res.success:
    (a0,b0)=res.x
  else:
    raise RuntimeError(res.message)
  l0=[x for x in l if x[4]>=cutoff]
  n=len(l0)
  sl=[]
  nit=1000
  for it in range(nit):
    m=[l0[randrange(n)] for i in range(n)]
    res=minimize(err,xx,args=(m,cutoff),method="SLSQP",bounds=bounds,options={"maxiter":1000})
    if res.success:
      (a,b)=res.x
      sl.append(b)
    else:
      raise RuntimeError(res.message)
  sl.sort()
  conf=0.9917
  (bmin,bmax)=sl[int((1-conf)/2*nit)],sl[int((1+conf)/2*nit)]
  print("Cutoff %2d : %3d LTLAs : Extra transmissibility %3.0f%% (%.0f - %.0f%%)"%(cutoff,n,b0*100,bmin*100,bmax*100))
