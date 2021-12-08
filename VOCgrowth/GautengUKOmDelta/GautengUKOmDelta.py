# Sources (as at 2021-12-06):
# https://sacoronavirus.co.za/latest-vaccine-statistics/
# https://coronavirus.data.gov.uk/details/vaccinations?areaType=overview&areaName=United%20Kingdom
# https://www.gov.uk/government/publications/covid-19-vaccine-weekly-surveillance-reports
# http://sonorouschocolate.com/covid19/index.php?title=Early_Omicron_Growth_Estimate
#
# See https://docs.google.com/spreadsheets/d/1wvQdozRGTpCvU44IQiT8D5z6dQVBR-9xTHlZ5vhFzW0/edit#gid=812895892
# for a spreadsheet illustration of this program

from math import exp,log
from random import normalvariate,random

# Input is vaccine rate array, attack rate. VR[i] = proportion of population who have had exactly i doses
# Assume vaccination status is independent of whether infected (yes, I know)
def IndImStatus(VR,AR):
  return [VR[0]*(1-AR),VR[0]*AR,
          VR[1]*(1-AR),VR[1]*AR,
          VR[2]*(1-AR),VR[2]*AR,
          VR[3]]

def getdoublingtime(
    # Vaccination rates
    # 0 dose, 1 dose, 2 dose, 3 dose
    VR_Gauteng=[0.430, 0.364, 0.206, 0.],
    VR_UK=[0.238, 0.068, 0.392, 0.302],
    
    # Attack rates
    AR_Gauteng=0.95,
    AR_UK=0.45,
    
    # Generation times
    T_Del=5,
    T_Om=5,

    # Growth rates
    g_Del_Gauteng=-0.043,
    g_Om_Gauteng=0.232,
    g_Del_UK=0.005,

    # Immunity coefficients vs new symptomatic infection
    # Map from immunity status as defined by
    # [Naive, pre-infected only, 1 dose, pre-infected + 1 dose, 2 dose, pre-infected + 2 dose, 3 dose]
    # Can also specify Im_Om as a logistic shift from Im_Del
    Im_Del=[0, 0.75, 0.30, 0.80, 0.60, 0.92, 0.92],
    Im_Om=[0, 0.15, 0.05, 0.40, 0.30, 0.70, 0.80],
    
    # Whether to print stuff
    pr=False
  ):
  
  assert abs(sum(VR_Gauteng)-1)<1e-6 and abs(sum(VR_UK)-1)<1e-6
  ImStatus_Gauteng=IndImStatus(VR_Gauteng,AR_Gauteng)
  ImStatus_UK=IndImStatus(VR_UK,AR_UK)

  if pr:
    print('                0d      i     1d   1d+i     2d   2d+i     3d')
    print('vs Delta:  ',end='')
    for x in Im_Del:
      print(' %5.1f%%'%(x*100),end='')
    print()
    print('vs Omicron:',end='')
    for x in Im_Om:
      print(' %5.1f%%'%(x*100),end='')
    print()
  
  RR_Gauteng_Del=1-sum(x*y for (x,y) in zip(ImStatus_Gauteng,Im_Del))
  RR_Gauteng_Om=1-sum(x*y for (x,y) in zip(ImStatus_Gauteng,Im_Om))
  RR_UK_Del=1-sum(x*y for (x,y) in zip(ImStatus_UK,Im_Del))
  RR_UK_Om=1-sum(x*y for (x,y) in zip(ImStatus_UK,Im_Om))
  
  RRR_Gauteng=RR_Gauteng_Om/RR_Gauteng_Del
  Rrel_obs_Gauteng=exp(T_Om*g_Om_Gauteng-T_Del*g_Del_Gauteng)
  RR0=Rrel_obs_Gauteng/RRR_Gauteng
  RRR_UK=RR_UK_Om/RR_UK_Del
  Rrel_pred_UK=RR0*RRR_UK
  # Rrel_pred_UK=exp(T_Om*g_Om_UK-T_Del*g_Del_UK)
  g_Om_UK=(log(Rrel_pred_UK)+T_Del*g_Del_UK)/T_Om
  return log(2)/g_Om_UK,RR0

def logistmult(p,x):
  if type(p)==list: return [logistmult(y,x) for y in p]
  return x*p/(x*p+1-p)

def logistadd(p,x):
  return logistmult(p,exp(x))

# Choose a random probability given the 95% confidence interval is [p0,p1]. (Assume the log odds are normal.)
def getrandprob(p0,p1):
  if p1==0: return 0
  if p0==1: return 1
  zconf=1.96
  x0=log(p0/(1-p0))
  y0=log(p1/(1-p1))
  mu=(x0+y0)/2
  sd=(y0-x0)/2/zconf
  z=normalvariate(mu,sd)
  return 1/(1+exp(-z))

def getrandefficacies():
  i0=getrandprob(.72,.87)
  d1=getrandprob(.52,.64)
  d2=getrandprob(.71,.78)
  d3=getrandprob(.92,.96)
  i1=logistadd(max(i0,d1),0.5*random())
  i2=logistadd(max(i1,d2),0.5+0.5*random())
  Im_Del=[0,i0,d1,i1,d2,i2,d3]
  Im_Om=logistadd(Im_Del,-2.3)
  return Im_Del,Im_Om
  
#Im_Del=[0, 0.75, 0.30, 0.80, 0.60, 0.92, 0.92]
#Im_Om=logistadd(Im_Del,normalvariate(-2.5,0.3))
#print(getdoublingtime(Im_Del=Im_Del,Im_Om=Im_Om))

def confint(l,conf):
  n=len(l)
  return (l[n//2],l[int(n*(1-conf)/2)],l[int(n*(1+conf)/2)])

nit=10000
pr=0
histdoub=[]
histRR0=[]
for i in range(nit):
  Im_Del,Im_Om=getrandefficacies()
  AR_Gauteng=getrandprob(0.85,0.98)
  AR_UK=getrandprob(0.40,0.50)
  doub,RR0=getdoublingtime(AR_Gauteng=AR_Gauteng,AR_UK=AR_UK,Im_Del=Im_Del,Im_Om=Im_Om)
  if pr:
    for p in Im_Del: print('   %5.1f'%(p*100),end='')
    print()
    for p in Im_Om: print('   %5.1f'%(p*100),end='')
    print()
    print(doub,RR0)
    print()
  histdoub.append(doub)
  histRR0.append(RR0)
histdoub.sort()
histRR0.sort()
print("Doubling period: %.2f (%.2f - %.2f) days"%confint(histdoub,0.95))
print("R0 ratio: %.2f (%.2f - %.2f)"%confint(histRR0,0.95))
