from stuff import *
import numpy as np
from math import  log

mindate='2021-05-15'

ladcsv=loadcsv('Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv')
sanger=loadcsv('lineages_by_ltla_and_week.tsv',sep='\t')

lad2reg=dict(zip(ladcsv['LAD19CD'],ladcsv['RGN19NM']))
regions=sorted(list(set(ladcsv['RGN19NM'])))
nreg=len(regions)

# Alpha, Delta maps from region to date to count
A={r:{} for r in regions}
D={r:{} for r in regions}
dates=set()
for (date,lad,lin,n) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
  if date>=mindate:
    reg=lad2reg[lad]
    if lin=='B.1.1.7': A[reg][date]=A[reg].get(date,0)+n
    if lin=='B.1.617.2': D[reg][date]=D[reg].get(date,0)+n
    dates.add(date)

for reg in regions:
  W=X=Y=XX=XY=YY=0
  for date in dates:
    x=datetoday(date)-datetoday(mindate)
    a,d=A[reg].get(date,0),D[reg].get(date,0)
    if a>0 and d>0:
      y=log(d/a)
      w=1/(1/a+1/d)
      W+=w
      X+=w*x
      Y+=w*y
      XX+=w*x*x
      XY+=w*x*y
      YY+=w*y*y
  M=np.array([[W, X], [X, XX]])
  R=np.array([Y,XY])
  S=np.linalg.solve(M,R)
  print("%-30s"%reg,-S[0]/S[1],S[1])
