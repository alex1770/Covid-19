from stuff import *
import sys
from math import sqrt,log
import numpy as np
from scipy.optimize import minimize

# Get SGTF/S-gene from fig. 13 from this
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/986378/Variants_of_Concern_Technical_Briefing_11_Data_England__1_.xlsx
# whose containing page is
# https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201

cases=loadcsv("ltla_2021-05-15.csv")
sgtf=loadcsv("Variants_of_Concern_Technical_Briefing_11_Data_England-1-fig13.csv")

date0='2021-04-25'
date1='2021-05-02'
date2='2021-05-09'
minsgeneclass=30

cases0={}
cases1={}
for (ac,date,n) in zip(cases['areaCode'],cases['date'],cases['newCasesBySpecimenDate']):
  if date>=date0 and date<date1: cases0[ac]=cases0.get(ac,0)+n
  if date>=date1 and date<date2: cases1[ac]=cases1.get(ac,0)+n

ltlas=set()
sgtfnum={}
for (x,n,r) in zip(sgtf['ltla.code'],sgtf['Number classifiable cases'],sgtf['Number S gene cases']):
  if n>=minsgeneclass and x in cases0 and x in cases1:
    ltlas.add(x)
    sgtfnum[x]=(n,r)
print("Using",len(ltlas),"LTLAs",file=sys.stderr)
ltlas=sorted(list(ltlas))

def LL(Q,R,X,Y):
  tot=0
  blah=0
  for ltla in ltlas:
    a=cases0[ltla]+X
    b=cases1[ltla]
    (c,d)=sgtfnum[ltla]
    A=a*(R-Q)+c
    x=Q/(R-Q)
    B=A*(x+Y)-b-d
    C=A*x*Y-b*Y-x*d
    s=sqrt(B**2-4*A*C)
    p=(-B+s)/(2*A)
    #if p<0 or p>1: print(ltla,p,a,b,c,d)
    p=min(max(p,0),1)
    tot+=-(a*(R-Q)+c)*p-c*Y-a*Q+b*log(a*(p*R+(1-p)*Q))+(d*log(p+Y) if d>0 else 0)
    blah+=b*log(cases0[ltla])
  #print(blah)
  return tot-blah

def MLL(xx): return -LL(*xx)

xx=[1,3,0,0]
bounds=[(0.1,1.4),(1.6,10),(0,100),(0,1)]
res=minimize(MLL,xx,method="SLSQP",bounds=bounds,options={"maxiter":1000})

print(res)
(Q,R,X,Y)=res.x

for ltla in ltlas:
  a=cases0[ltla]+X
  b=cases1[ltla]
  (c,d)=sgtfnum[ltla]
  A=a*(R-Q)+c
  x=Q/(R-Q)
  B=A*(x+Y)-b-d
  C=A*x*Y-b*Y-x*d
  s=sqrt(B**2-4*A*C)
  p=(-B+s)/(2*A)
  #if p<0 or p>1: print(ltla,p,a,b,c,d)
  p=min(max(p,0),1)
  def f(p):
    return -(a*(R-Q)+c)*p-c*Y-a*Q+b*log(a*(p*R+(1-p)*Q))+(d*log(p+Y) if d>0 else 0)
  break
  #blah+=b*log(cases0[ltla])
  print(ltla,"%5.3f %5.3f   %5.3f   %5.1f %5d %5d %5d"%(p+Y,d/c,p*R+(1-p)*Q,a,b,c,d))
  #print(d/c,b/a)
