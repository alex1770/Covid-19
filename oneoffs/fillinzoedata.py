#/usr/bin/env python3

import os,ast,json,sys
from os.path import join
from collections import defaultdict

floatkeys=["long", "lat", "st_areasha", "st_lengths", "corrected_covid_positive", "cases", "cases_pm", "percentage", "discrete_percentage"]
intkeys=["cartodb_id", "objectid", "bng_e", "bng_n", "respondent", "predicted_covid_positive_count", "population", "discrete_cases_pm"]

tdir="zoemapdata"
l=sorted((x for x in os.listdir(tdir) if x[:2]=='20'), reverse=True)

#outdir=tdir
outdir="tempd1"

# Categories
cats=["respondent", "population", "corrected_covid_positive"]

ee={}
lookup=defaultdict(dict)
for date in l:
  with open(join(tdir,date),'r') as fp:
    dd=json.load(fp)
  # dd :: location -> category -> value
  # ee :: date -> category -> location -> value
  # lookup :: category -> location -> date
  ee[date]={}
  for c in cats:
    ee[date][c]={}
    for loc in dd:
      if c in dd[loc]:
        ee[date][c][loc]=dd[loc][c]
        lookup[c][loc]=date
  for c in cats:
    for loc in dd:
      if c not in dd[loc]:
        #print("Missing",date,loc,c)
        if loc not in lookup[c]:
          print(date,loc,c,"not fillinable from future")
        else:
          date1=lookup[c][loc]
          print("Filling in",date,loc,c,"from",date1)
          meet=set(ee[date][c]).intersection(ee[date1][c])
          tot=sum(ee[date][c][m] for m in meet)
          tot1=sum(ee[date1][c][m] for m in meet)
          dd[loc][c]=tot/tot1*ee[date1][c][loc]
          if c in intkeys: dd[loc][c]=int(dd[loc][c]+0.5)
  with open(join(outdir,date),'w') as fp:
    json.dump(dd,fp,indent=2)
