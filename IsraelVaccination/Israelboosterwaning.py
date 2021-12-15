from stuff import *
from random import sample
from math import sqrt

# Whether to do sensitivity analysis
sensitivity=False

# From https://datadashboard.health.gov.il/COVID-19/general
# https://data.gov.il/dataset/covid-19/resource/57410611-936c-49a6-ac3c-838171055b1f
# https://data.gov.il/dataset/covid-19/resource/9b623a64-f7df-4d0c-9f57-09bd99a88880
vax=loadcsv('vaccinated-per-day-2021-11-23.csv')
cases=loadcsv('Israel-cases-among-vaccinated-164.csv')
ages=sorted(list(set(vax['age_group'])))
assert ages==sorted(list(set(cases['Age_group'])))
  
# Population (thousands) by 5-year age bands from WPP2019's 2020 estimate
# https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx
WPPpop=[848, 823, 738, 667, 624, 573, 565, 558, 544, 497, 418, 378, 350, 354, 289, 171, 124, 84, 40, 12, 1]
totpopwpp={}
for age in ages:
  if age[-1]=='+':
    a=int(age[:-1])
    v=sum(WPPpop[a//5:])
  else:
    a,b=map(int,age.split('-'))
    v=sum(WPPpop[a//5:(b+1)//5])
  totpopwpp[age]=v*1000
#  0-19  3076000
# 20-29  1197000
# 30-39  1123000
# 40-49  1041000
# 50-59   796000
# 60-69   704000
# 70-79   460000
# 80-89   208000
#   90+    53000

# Sarah's population estimate = max_{d=1,2,3} (d-dose vaccinated)/(d-dose proportion vaccinated)   (+tiny adjustment to make percentages add to 100%)
totpopsarah={'20-29':1321161, '30-39':1214249, '40-49':1114305, '50-59':870068, '60-69':747070, '70-79':512664, '80-89':233857, '90+':55928}

#totpop=totpopsarah
totpop={age: sqrt(totpopwpp[age]*totpopsarah[age]) for age in totpopsarah}
al=0.6
totpop={age: totpopwpp[age]**al*totpopsarah[age]**(1-al) for age in totpopsarah}
  
print("Using populations:")
for age in totpop:
  print("%5s  %7d"%(age,totpop[age]))
print()

def getint(s):
  if type(s)==int: return s
  if s=='': return 0
  if s[0]=='<': return 1
  return int(s)

# totvax[ageband,date,dose] = total number vaccinated up to 'date'
totvax={}
totvaxr={(age,d):0 for age in ages for d in [1,2,3]}
for (date,age,d1,d2,d3) in zip(vax['VaccinationDate'], vax['age_group'], vax['first_dose'], vax['second_dose'], vax['third_dose']):
  totvaxr[age,1]+=getint(d1);totvax[age,date,1]=totvaxr[age,1]
  totvaxr[age,2]+=getint(d2);totvax[age,date,2]=totvaxr[age,2]
  totvaxr[age,3]+=getint(d3);totvax[age,date,3]=totvaxr[age,3]

weeks=sorted(week for week in set(cases['Week']) if week>='2021-09-19')

def tau(perm):
  n=len(perm)
  tau=0
  for i in range(n-1):
    for j in range(i+1,n):
      tau+=(perm[i]<perm[j])-(perm[i]>perm[j])
  return tau

def pvalue(perms,N=10000):
  t0=sum(tau(perm) for perm in perms)
  nn=[len(perm) for perm in perms]
  m=0
  for i in range(N):
    qerms=[sample(range(n),n) for n in nn]
    t=sum(tau(sample(range(n),n)) for n in nn)
    m+=(t>=t0)
  return m/N

offset=20
infinity=datetoday('2021-11-23')
if sensitivity:
  popmults=[0.95, 1, 1.05]
  weekdayoffsets=[0,3,6,'average','infinity']
else:
  popmults=[1]
  weekdayoffsets=['average']
for popmult in popmults:
  for weekdayoffset in weekdayoffsets:
    RRLL=[]
    for age in ages[1:-1]:
      print('Age band:',age)
      print('Weekdayoffset:',weekdayoffset)
      print('Population multiplier: %g'%popmult)
      print('VC=vaccinated cases, VP=vaccinated population, NVC=non-vaccinated cases, NVP=non-vaccinated population, RR=(VC/VP)/(NVC/NVP)')
      print("Week of cases               VC      VP      NVC     NVP       RR    1-RR    VC/NVC")
      RRL=[]
      for (week,a,vaxcases,nonvaxcases) in zip(cases['Week'], cases['Age_group'], cases['positive_above_20_days_after_3rd_dose'], cases['Sum_positive_without_vaccination']):
        if a!=age or week not in weeks: continue
        vaxcases=getint(vaxcases)
        nonvaxcases=getint(nonvaxcases)
        if weekdayoffset=='infinity':
          days=[infinity]
        elif weekdayoffset=='average':
          d=datetoday(week[:10])
          days=range(d,d+7)
        else:
          d=datetoday(week[:10])
          days=[d+weekdayoffset]
        vaxpop1=sum(totvax[age,daytodate(d),1] for d in days)/len(days)
        vaxpop3=sum(totvax[age,daytodate(d-offset),3] for d in days)/len(days)
        nonvaxpop=totpop[age]*popmult-vaxpop1
        RR=vaxcases/max(vaxpop3,1e-10)/(nonvaxcases/nonvaxpop)
        RRL.append(RR)
        if vaxpop3>0:
          print(week," %5d %7d    %5d %7d   %5.1f%%  %5.1f%%    %6.4f"%(vaxcases,vaxpop3,nonvaxcases,nonvaxpop,RR*100,(1-RR)*100,vaxcases/nonvaxcases))
      print("Age=%s, weekdayoffset=%s, popmult=%g, Tau test p-value %.3g"%(age,str(weekdayoffset),popmult,pvalue([RRL],20000)))
      print()
      RRLL.append(RRL)
    print("Overall tau test p-value %.5f at popmult=%g, weekdayoffset"%(pvalue(RRLL,100000),popmult),weekdayoffset)
    print()
  
