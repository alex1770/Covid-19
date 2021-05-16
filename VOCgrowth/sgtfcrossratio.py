from stuff import *
from scipy.stats import gamma

# Get SGTF/S-gene from fig. 15 from this
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/986378/Variants_of_Concern_Technical_Briefing_11_Data_England__1_.xlsx
# whose containing page is
# https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201

sgtftab=loadcsv("Variants_of_Concern_Technical_Briefing_11_Data_England__1_fig15.csv")

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

mgt=5# Mean generation time
for weekpair in [('21/04/2021','28/04/2021',7),('28/04/2021','05/05/2021',7),('21/04/2021','05/05/2021',14)]:
  print("Comparing weeks (w/c)",weekpair[0],"and",weekpair[1])
  # Get counts into the form [1st week SGTF, 1st week S gene, 2nd week SGTF, 2nd week S gene]
  sgtf={x:[0,0,0,0] for x in sgtftab['Region']+['England total']}
  for (r,w,s,c) in zip(sgtftab['Region'], sgtftab['week'], sgtftab['S gene detection result'], sgtftab['case count']):
    if w in weekpair[:2]:
      i=2*(w==weekpair[1])+(s=='S gene positive')
      sgtf[r][i]=c
      sgtf['England total'][i]+=c

  print("                           ----Cross ratio----        -----R factor------")
  print("                           Low     Med    High        Low     Med    High")
  regions=sorted(list(set(sgtftab['Region'])))+['England total']
  days=weekpair[2]
  for r in regions:
    [a,b,c,d]=sgtf[r]
    (low,med,high)=credint(a,b,c,d,0.98)
    f=mgt/days
    print("%-22s  %6.2f  %6.2f  %6.2f     %6.2f  %6.2f  %6.2f"%(r,low,med,high,low**f,med**f,high**f))
    
  print()
