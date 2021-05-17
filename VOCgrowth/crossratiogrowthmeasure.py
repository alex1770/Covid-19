from stuff import *
from scipy.stats import gamma,norm
from math import exp,sqrt

source="COG"
#source="Sanger"

print("Using data source:",source)
if source=="Sanger":
  # Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv
  sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
elif source=="COG":
  # https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv
  cog=loadcsv("cog_metadata.csv")
else:
  raise RuntimeError("Unrecognised source: "+source)
  
regiondat=loadcsv("Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.csv")

VV=["B.1.1.7","B.1.617.2"]# Variants to compare

mgt=5# Mean generation time

mindate='2021-04-10'

def credint(a,b,c,d,conf=0.95):
  nit=1000000
  eps=0.5
  A=gamma.rvs(a+eps,size=nit)
  B=gamma.rvs(b+eps,size=nit)
  C=gamma.rvs(c+eps,size=nit)
  D=gamma.rvs(d+eps,size=nit)
  l=A*D/(B*C)
  l.sort()
  return (l[int((1-conf)/2*nit)], l[nit//2], l[int((1+conf)/2*nit)])

def confint(a,b,c,d,conf=0.95):
  OR=(a*d)/(b*c)
  z=norm.ppf((1+conf)/2)
  m=exp(z*sqrt(1/a+1/b+1/c+1/d))
  return [OR/m,OR,OR*m]

print("Comparing",VV[0],"and",VV[1])

ltla2region={}
for (ltla,reg) in zip(regiondat['LAD19CD'],regiondat['RGN19NM']):
  ltla2region[ltla]=reg

if source=="Sanger":
  regions=sorted(list(set(ltla2region.values())))+['England']
  dates=sorted([date for date in set(sanger['WeekEndDate']) if date>=mindate])
  data={region:{} for region in regions}
  for region in regions:
    data[region]={}
    for date in dates:
      data[region][date]={v:0 for v in VV}
  for (date,ltla,lineage,count) in zip(sanger['WeekEndDate'],sanger['LTLA'],sanger['Lineage'],sanger['Count']):
    if date>=mindate and lineage in VV:
      for region in [ltla2region[ltla],'England']:
        data[region][date][lineage]=data[region][date].get(lineage,0)+count
elif source=="COG":
  regions=sorted(list(set(cog['adm1'])))+['UK']
  weekenddate='2021-05-08'# Use this to align to Sanger weeks. Weeks end on Saturday
  weekendday=datetoday(weekenddate)
  minday=datetoday(mindate)
  minday+=(weekendday+1-minday)%7
  maxday=datetoday(max(list(set(cog['sample_date']))))
  maxday-=(maxday-weekendday)%7
  dates=[daytodate(day) for day in range(minday-1,maxday+7,7)]
  data={region:{} for region in regions}
  for region in regions:
    data[region]={}
    for date in dates:
      data[region][date]={v:0 for v in VV}
  for (country,adm1,date,lineage) in zip(cog['country'],cog['adm1'],cog['sample_date'],cog['lineage']):
    if country=="UK" and lineage in VV:
      day=datetoday(date)
      if day>=minday and day<=maxday:
        date=daytodate(day+(maxday-day)%7)
        for region in [adm1,"UK"]:
          data[region][date][lineage]=data[region][date].get(lineage,0)+1
  
for region in regions:
  print(region)
  dat=data[region]
  C=D=0
  for i in range(len(dates)):
    d=dat[dates[i]]
    A,B=C,D
    C,D=d[VV[0]],d[VV[1]]
    print(dates[i],"%6d %6d"%(C,D),end=' ')
    if A==0 or B==0 or C==0 or D==0: print();continue
    (low,med,high)=[x**(mgt/7) for x in confint(A,B,C,D)]
    print("%6.1f%% ( %6.1f - %6.1f%% )"%((med-1)*100,(low-1)*100,(high-1)*100))
  print()
