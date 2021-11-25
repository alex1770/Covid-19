from stuff import *
from random import sample

# From https://datadashboard.health.gov.il/COVID-19/general
# https://data.gov.il/dataset/covid-19/resource/57410611-936c-49a6-ac3c-838171055b1f
# https://data.gov.il/dataset/covid-19/resource/9b623a64-f7df-4d0c-9f57-09bd99a88880
vax=loadcsv('vaccinated-per-day-2021-11-23.csv')
cases=loadcsv('Israel-cases-among-vaccinated-164.csv')
ages=sorted(list(set(vax['age_group'])))
assert ages==sorted(list(set(cases['Age_group'])))
  
if 0:
  # Population (thousands) by 5-year age bands from WPP2019's 2020 estimate
  # https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx
  WPPpop=[848, 823, 738, 667, 624, 573, 565, 558, 544, 497, 418, 378, 350, 354, 289, 171, 124, 84, 40, 12, 1]
  
  totpop={}
  print("Using populations")
  for age in ages:
    if age[-1]=='+':
      a=int(age[:-1])
      v=sum(WPPpop[a//5:])
    else:
      a,b=map(int,age.split('-'))
      v=sum(WPPpop[a//5:(b+1)//5])
    totpop[age]=v*1000*1.02

# Sarah's estimate = max_{d=1,2,3} (d-dose vaccinated)/(d-dose proportion vaccinated)   (+tiny adjustment to make percentages add to 100%)
totpop={'20-29':1321161, '30-39':1214249, '40-49':1114305, '50-59':870068, '60-69':747070, '70-79':512664, '80-89':233857, '90+':55928}
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
  prev=daytodate(datetoday(date)-1)
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

def pvalue(perm,N=1000):
  t=tau(perm)
  n=len(perm)
  m=nit=0
  while m<N or nit-m<N:
    qerm=sample(range(n),n)
    m+=(tau(qerm)>=t)
    nit+=1
  return m/nit

for age in ages[1:-2]:
  print('Age band',age)
  print("Week of cases               VC      VP     NVC     NVP       RR    1-RR    VC/NVC")
  RRL=[]
  for (week,a,vaxcases,nonvaxcases) in zip(cases['Week'], cases['Age_group'], cases['positive_above_20_days_after_3rd_dose'], cases['Sum_positive_without_vaccination']):
    if a!=age or week not in weeks: continue
    vaxcases=getint(vaxcases)
    nonvaxcases=getint(nonvaxcases)
    prev=daytodate(datetoday(week[:10])-17)
    vaxpop=totvax[age,prev,3]
    nonvaxpop=totpop[age]-totvax[age,prev,1]
    RR=vaxcases/vaxpop/(nonvaxcases/nonvaxpop)
    RRL.append(RR)
    print(week,"  %4d %7d    %4d %7d   %5.1f%%  %5.1f%%    %6.4f"%(vaxcases,vaxpop,nonvaxcases,nonvaxpop,RR*100,(1-RR)*100,vaxcases/nonvaxcases))
  print("p-value %.3g"%(pvalue(RRL)))
  print()
