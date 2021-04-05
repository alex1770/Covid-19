import sys,csv,time,calendar
from math import log,sqrt
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln

keephalfterm=True
#maxdate='9999-99-99'
maxdate='2021-03-08'
agerange=(10,19)

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

# wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesLFDConfirmedPCRBySpecimenDate&metric=newCasesLFDOnlyBySpecimenDate&metric=newLFDTests&metric=newCasesBySpecimenDate&format=csv' -O engcasesbyspecimen.csv
# wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesByPublishDate&format=csv' -O casesbypublication.csv
# wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaName=England&metric=maleCases&metric=femaleCases&format=csv' -O casesbyage.csv

# 2020 UK population estimate by 5-year age band (in thousands)
# from https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2019_POP_F15_1_ANNUAL_POPULATION_BY_AGE_BOTH_SEXES.xlsx
pop=[3924, 4120, 3956, 3686, 4075, 4484, 4707, 4588, 4308, 4296, 4635, 4539, 3905, 3382, 3388, 2442, 1737, 1078, 491, 130, 16]

# Estimate England population by multiplying UK population by multiplying by 55.98/66.65
population=0
for (i,n) in enumerate(pop):
  if agerange[0]<=i*5 and agerange[1]+1>=(i+1)*5: population+=55.98/66.65*1000*n

def loadcsv(fn,keephalfterm=True,maxdate='9999-99-99'):
  dd={}
  with open(fn,"r") as fp:
    reader=csv.reader(fp)
    headings=[x.strip() for x in next(reader)]
    for row in reader:
      if (keephalfterm or row[0]!='2021-02-17') and row[0]<=maxdate:
        for (name,x) in zip(headings,row):
          x=x.strip()
          if x.isdigit(): x=int(x)
          dd.setdefault(name,[]).append(x)
  return dd

print("Using LFD school numbers from table 7 of \"Tests conducted\" spreadsheet, available from https://www.gov.uk/government/collections/nhs-test-and-trace-statistics-england-weekly-reports")
#https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/973128/tests_conducted_2021_03_25.ods")
if keephalfterm: print("Keeping half term 14-20 February")
else: print("Discarding half term 14-20 February")
print("Using school LFD data up to w/e",maxdate)
print("Using age range %d-%d, population %.1fm"%(agerange+(population/1e6,)))
dd=loadcsv("LFDschooltests.csv",keephalfterm=keephalfterm,maxdate=maxdate)
dd['LFDnum']=list(map(sum,zip(dd['LFDpos'],dd['LFDneg'])))# Denominator = positive tests + negative tests (ignore unknown/void tests)

# cc=loadcsv("engcasesbyspecimen.csv")
# datetocases=dict(zip(cc['date'],cc['newCasesBySpecimenDate']))
# print("Using England cases by specimen date")

# cc=loadcsv("casesbypublication.csv")
# datetocases=dict(zip(cc['date'],cc['newCasesByPublishDate']))
# print("Using UK cases by publication date from https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=newCasesByPublishDate&format=csv")
# print()

cc=loadcsv("casesbyage.csv")
cumdatetocases={}
for (date,metric,age,value) in zip(cc['date'],cc['metric'],cc['age'],cc['value']):
  if age[-1]=='+': ar=[int(age[:-1]),149]
  else: ar=[int(x) for x in age.split('_to_')]
  if ar[0]>=agerange[0] and ar[1]<=agerange[1]: cumdatetocases[date]=cumdatetocases.get(date,0)+value
print("Using England cases by age from https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaName=England&metric=maleCases&metric=femaleCases&format=csv")
print()

# Work out constant part of negative binomial denominator. It doesn't affect the optimisation, but it makes the log likelihood output more meaningful
LL0=sum(gammaln(a+1) for a in dd['LFDpos'])

# SLSQP seems to be happier if the variables we are optimising over are of order 1, so use scale factors
scale0=1e4
scale2=1e3

# Return -log(likelihood) assuming negative binomial distribution with common dispersion parameter r=xx[2]*scale2
def err(xx):
  r=xx[2]*scale2
  n=len(cases)
  LL=-n*gammaln(r)
  for (a,b,c) in zip(dd['LFDpos'],dd['LFDnum'],cases):
    mu=(xx[0]/scale0+xx[1]*c/population)*b
    p=mu/(mu+r)
    LL+=gammaln(r+a)+r*log(1-p)+a*log(p)
  return LL0-LL

# Return negative second derivative of log(likelihood) wrt xx[0]/scale0 (aka observed Fisher information)
def LL2(xx):
  eps=1e-3
  e1=err(xx)
  while 1:
    e2=err([xx[0]+eps*scale0]+list(xx[1:]))
    if abs(e2-e1)<1e-3*(e1+e2): break
    eps/=2
  e0=err([xx[0]-eps*scale0]+list(xx[1:]))
  return (e0-2*e1+e2)/eps**2


best=(-1e9,)
for offset in range(0,1):# Seems more principled to force offset 0 (what you'd expect it to be) rather than allow it as a free parameter
  cases=[]
  for date in dd["WeekEnding"]:
    day=datetoday(date)
    #cases.append(sum(datetocases[daytodate(d)] for d in range(day+offset-6,day+offset+1)))
    cases.append(cumdatetocases[daytodate(day+offset)]-cumdatetocases[daytodate(day+offset-7)])
  #cases=cases[1:]+[cases[-1]*0.95]#alter
  if 0:# Sensitivity check, to check that the correlation between LFDpos and cases is not due to LFDpos feeding directly into cases, or possibly a correction term
    for i in range(len(cases)): cases[i]-=dd['LFDpos'][i]*0.3
  
  res=minimize(err,[1,1,1],method="SLSQP",bounds=[(1e-9,scale0),(1e-9,100),(0.01,100)],options={"maxiter":1000})
  if not res.success: raise RuntimeError(res.message)
  LL=-res.fun
  xx=res.x
  fisher=LL2(xx)
  print("Offset %3d. Log likelihood = %8.3f. al = %6.4f%% . be = %8.3f . r = %8.3g"%((offset,LL,xx[0]/scale0*100,xx[1],xx[2]*scale2)))
  if 1:
    for (a,b,c,w) in zip(dd['LFDpos'],dd['LFDnum'],cases,dd['WeekEnding']):
      lam=(xx[0]/scale0+xx[1]*c/population)*b
      print("   w/e %s %8.3f  %6d"%(w,a*log(lam)-lam-gammaln(a+1),c))
    print()
  if LL>best[0]: best=(LL,offset,res,fisher,cases)
print()

(LL,offset,res,fisher,cases)=best
print("Best offset %d. Log likelihood = %g"%(offset,LL))
date=dd["WeekEnding"][-1]
day=datetoday(date)-6
print("Offset %d means LFD numbers for week %s - %s are related to national case numbers in the week %s - %s"%(offset,daytodate(day),date,daytodate(day+offset),daytodate(day+offset+6)))
xx=res.x
fpr=xx[0]/scale0*100
error=1.96/sqrt(fisher)*100
print("Best estimate: LFDpos/LFDnum = %.3g + %.3g*(weekly case rate)"%(xx[0]/scale0,xx[1]))
print("               where (weekly case rate) = (number of confirmed cases in a week) / %g"%population)
print("False positive rate estimate: %.2g%% (%.2g%% - %.2g%%)"%(fpr,fpr-error,fpr+error))

data=[]
with open("graph","w") as fp:
  r=xx[2]*scale2
  for (numlfdpos,numlfdtests,ncases) in zip(dd['LFDpos'],dd['LFDnum'],cases):
    sd=1.96*sqrt(numlfdpos*(numlfdpos/r+1))
    data.append((numlfdpos,numlfdtests,ncases,sd))
    print("%6d %8d %6d %9.3f"%(data[-1]),file=fp)

from subprocess import Popen,PIPE
po=Popen("gnuplot",shell=True,stdin=PIPE)
pp=po.stdin

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): pp.write((' '.join(map(str,s))+'\n').encode('utf-8'))

write('set terminal pngcairo font "sans,13" size 1920,1280')
write('set bmargin 5;set lmargin 15;set rmargin 10;set tmargin 5')
write('set title "Comparison of LFD positivity rate among secondary schools students with the rate of age %d-%d confirmed cases across England"'%agerange)
write('set xtics nomirror')
write('set xlabel "Total new confirmed cases age %d-%d in England over week-long period, as a percentage of the age %d-%d population"'%(agerange+agerange))
write('set ylabel "Total new LFD positive tests in secondary school students as a percentage of tests taken"')

outfn="lfdcases.png"
write('set output "%s"'%outfn)
write('plot [0:] [0:] "-" using 1:2 lw 12 title "%%LFD student positivity over a week vs %%new cases age %d-%d in England in same week"'%agerange)
for (numlfdpos,numlfdtests,ncases,sd) in data:
  write(100*ncases/population,100*numlfdpos/numlfdtests)
write('e')
print("Written graph to %s"%outfn)

outfn="lfdcaseswithfit.png"
write('set output "%s"'%outfn)
write('plot [0:] [0:] "-" using 1:2:3 with errorbars lw 3 title "%%LFD student positivity over a week vs %%new cases age %d-%d in England in same week", %g+%g*x lw 3'%(agerange+(100*xx[0]/scale0,xx[1])))
r=xx[2]*scale2
for (numlfdpos,numlfdtests,ncases,sd) in data:
  sd=1.96*sqrt(numlfdpos*(numlfdpos/r+1))
  write(100*ncases/population,100*numlfdpos/numlfdtests,100*sd/numlfdtests)
write('e')
pp.close();po.wait()
print("Written graph to %s"%outfn)

