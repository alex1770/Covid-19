# Feasibility check to see how possible it is to reverse-engineer individual week's
# figures from 4-week averages in vaccine surveillance report.

from stuff import *

pmd=loadcsv("vaccine-surveillance-reports/vaccine surveillance.csv")
ages=sorted(list(set(pmd['Age'])))
types=sorted(list(set(pmd['Type'])))
outcomes=['Total ', 'Unlinked', 'Not vaccinated', 'Received one dose 1-20 days before specimen date', 'Received one dose, 21 days or more before specimen date', 'Second dose 14 days or more before specimen date']

for age in ages:
  for ty in types:
    for outcome in outcomes:
      series=[v for (a,t,v) in zip(pmd['Age'],pmd['Type'],pmd[outcome]) if a==age and t==ty]

      n=len(series)
      
      lbs=[]
      for i in range(3):
        lb=0
        maxlb=0
        for j in range(2+i,n+1,4):
          lb+=series[j-2]-series[j-1]
          maxlb=max(maxlb,lb)
        lbs.append(maxlb)
        #print("x%d >="%i,maxlb)
        
      ub=minub=series[0]
      for j in range(5,n+1,4):
        ub+=series[j-1]-series[j-2]
        minub=min(minub,ub)
      #print("x0+x1+x2 <=",minub)

      if sum(lbs)>minub: print("Impossible",end='')
      elif minub>0: print("%10.3f"%(1-sum(lbs)/minub),end='')
      else: print("%10.3f"%0,end='')
      print("  %10s  %15s  "%(age,ty),outcome,'  ',series)
