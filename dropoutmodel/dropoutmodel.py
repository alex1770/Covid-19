# Starting from ONS dropout data, try to make model predicting proportions of dropouts of OR, N, S.
# Following on from Theo Sanderson's analysis at https://theo.io/post/2021-01-22-ons-data/.
# Data from tabs 6a, 6b of spreadsheet here https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/bulletins/coronaviruscovid19infectionsurveypilot/22january2021/relateddata
# as transcribed at https://github.com/theosanderson/theo.io/blob/master/content/post/2021-01-22-ons-data/ons_ct.csv

# Model:
# Let logistic(x)=exp(x)/(1+exp(x))
# Prevalence of B.1.1.7 = logistic(g*(t-l(Region)))
# Choose a uniform random number Z in [-5,5], which is fixed for the three genes,
# then probability of dropout for gene X (N, OR or S) = logistic(b*(Ct+Z-a_X)).
# Parameters (14):
#   a_X         : 3 parameters, one for each gene N, OR and S, encoding their "robustness" (lower = more fragile)
#   b           : 1 parameter encoding dependence of dropout probability on Ct
#   g           : 1 parameter for the relative growth of B.1.1.7 compared to other variants
#   l(Region)   : 9 parameters, one for each region, encoding takeover time of B.1.1.7

import sys,time,calendar,csv
from math import log,exp,sqrt
import numpy as np
from scipy.optimize import minimize

date0="2020-09-01"

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

class testdata:
  def __init__(self,date):
    self.date=date
    self.t=datetoday(date)-datetoday(date0)
    self.p=np.zeros([2,2,2])# p[whether N][whether OR][whether S] = proportion (adding to 1)
    self.Ct=0
  def __repr__(self):
    s=self.date+" :"
    for t in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(1,1,1)]: s+=" %5.1f"%(self.p[t]*100)
    s+=" : %.1f"%(self.Ct)
    return s

# csv headings:
# 0             1       2               3       4       5       6       7       8       9       10
# RegionType	Region	Week started	N only	OR only	S only	OR+N	OR+S	N+S	OR+N+S	Mean
#
# 11               12                   13              14              15
# 10th Percentile  25th Percentile	50th Percentile	75th Percentile	90th Percentile

data={}
with open("ons_ct.csv","r") as fp:
  reader=csv.reader(fp)
  headings=next(reader)
  for row in reader:
    if row[0]=="EnglandRegion":
      d=testdata(time.strftime("%Y-%m-%d",time.strptime(row[2],"%d %B %Y")))
      d.p[1][0][0]=float(row[3])# N
      d.p[0][1][0]=float(row[4])# OR
      d.p[0][0][1]=float(row[5])# S
      d.p[1][1][0]=float(row[6])# OR+N
      d.p[0][1][1]=float(row[7])# OR+S
      d.p[1][0][1]=float(row[8])# N+S
      d.p[1][1][1]=float(row[9])# OR+N+S
      d.p=d.p/d.p.sum()
      d.Ct=float(row[10])# Mean Ct
      data.setdefault(row[1],[]).append(d)

regions=sorted(list(data))

# Initial parameter values and bounds
xx=np.zeros(14);bounds=[[None,None] for i in range(14)]
xx[0]=xx[1]=xx[2]=30;  bounds[0]=bounds[1]=bounds[2]=(10,50)
xx[3]=1.0;             bounds[3]=(-1,3)
xx[4]=0.05;            bounds[4]=(-0.1,0.5)
for i in range(5,14): xx[i]=90;bounds[i]=(30,180)

# Work out expected dropout matrix
def estimatedropoutmatrix(r,d,xx):
  tc=np.zeros([2,2,2])
  for offset in range(-5,6):# Integrate over viral load
    p=1/(1+exp(-xx[4]*(d.t-xx[5+r])))# Relative prevalence of B.1.1.7
    dp=[1/(1+exp(-xx[3]*(d.Ct-a+offset*1.0))) for a in xx[:3]]# Probability of dropout for each gene
    dp[2]=p+(1-p)*dp[2]# Can treat B.1.1.7 as something that increases probability of S gene dropout
    c=np.zeros([2,2,2])
    for i in range(2):
      for j in range(2):
        for k in range(2):
          c[i,j,k]=(dp[0] if i==0 else 1-dp[0])*(dp[1] if j==0 else 1-dp[1])*(dp[2] if k==0 else 1-dp[2])
    tc+=c
  tc[0,0,0]=0;tc=tc/tc.sum()# Remove all-dropout and renormalise
  return tc
      
# Calculate modelling error associated with parameters xx[]
def err(xx):
  E=0
  for (r,region) in enumerate(regions):
    #if region!="London": continue
    for d in data[region]:
      c=estimatedropoutmatrix(r,d,xx)
      # Compare estimated and actual dropout matrices
      E+=((c-d.p)**2).sum()
  return E

res=minimize(err,xx,method="SLSQP",bounds=bounds,options={"maxiter":1000})

print(res.message)
if not res.success: sys.exit(1)
n=sum(len(data[region]) for region in regions)*7
print("RMS error = %.1f%%"%(sqrt(res.fun/n)*100))
print()

now=time.strftime('%Y-%m-%d',time.localtime())
tnow=datetoday(now)-datetoday(date0)
xx=res.x
print("Robustness of N  = %.1f"%xx[0])
print("Robustness of OR = %.1f"%xx[1])
print("Robustness of S  = %.1f"%xx[2])
print("Dependence of dropout on Ct = %.3f"%xx[3])
print("Relative growth rate per day of B.1.1.7 = %.2f%%"%(xx[4]*100))
print()
print("Region                    Est'd crossover     Extrapolated %relative")
print("------                    date of B.1.1.7     prevalence on",now)
print("                          ---------------     ------------------------")
for (r,region) in enumerate(regions):
  p=1/(1+exp(-xx[4]*(tnow-xx[5+r])))
  print("%-24s  %s        %6.1f"%(region,daytodate(datetoday(date0)+xx[5+r]+.5),p*100))
print()

def printdropouts(m):
  for t in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(1,1,1)]:
    print(" %5.1f"%(m[t]*100),end="")

print("                                                           Actual                                        Estimate                 ")
print("                                         ------------------------------------------     ------------------------------------------")
print("                  Region       Date         N    OR     S  OR+N  OR+S   N+S  OR+N+S        N    OR     S  OR+N  OR+S   N+S  OR+N+S")
for (r,region) in enumerate(regions):
  #if region!="London": continue
  for d in data[region]:
    print("%24s"%region,d.date,end="    ")
    printdropouts(d.p)
    print("     ",end="")
    c=estimatedropoutmatrix(r,d,xx)
    printdropouts(c)
    print()
