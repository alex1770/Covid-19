from stuff import *
from scipy.stats import gamma,norm
from math import exp,sqrt

source="COG"
#source="Sanger"
#source="PHE"
#source="SGTF"


print("Using data source:",source)
if source=="Sanger":
  # Sanger data from https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv
  sanger=loadcsv("lineages_by_ltla_and_week.tsv",sep='\t')
elif source=="COG":
  # https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv
  cog=loadcsv("cog_metadata.csv")
elif source=="PHE" or source=="SGTF":
  pass
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
  maxday=datetoday(max(list(set(cog['sample_date']))))-3
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
elif source=="PHE":
  regions=['England']
  # From VOC government data: https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data
  dat={
    '2021-04-28':{'B.1.1.7':8466 , 'B.1.617.2':202 },
    '2021-05-05':{'B.1.1.7':6795 , 'B.1.617.2':318 },
    '2021-05-12':{'B.1.1.7':9141 , 'B.1.617.2':793 },
    '2021-05-17':{'B.1.1.7':7066 , 'B.1.617.2':2111 }
  }
  data={'England':dat}
  dates=sorted(list(dat))
elif source=="SGTF":
  regions=['England']
  # From fig 14 of https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/986469/Variants_of_Concern_Technical_Briefing_11_Data_England-1.xlsx
  dat={
    '2021-03-24':{'B.1.1.7':14079, 'B.1.617.2':118 },
    '2021-03-31':{'B.1.1.7':9640 , 'B.1.617.2':186 },
    '2021-04-07':{'B.1.1.7':6923 , 'B.1.617.2':197 },
    '2021-04-14':{'B.1.1.7':5548 , 'B.1.617.2':324 },
    '2021-04-21':{'B.1.1.7':4558 , 'B.1.617.2':460 },
    '2021-04-28':{'B.1.1.7':3017 , 'B.1.617.2':866 },
    '2021-05-05':{'B.1.1.7':2642 , 'B.1.617.2':1632 },
  }
  data={'England':dat}
  dates=sorted(list(dat))



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
