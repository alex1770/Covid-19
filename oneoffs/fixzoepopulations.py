#/usr/bin/env python3

# Populations (used as a denominator for cases-per-million) were significantly incorrect
# (too low) for all regions before 2020-09-30. This program rewrites population values
# prior to 2020-09-30 to their 2020-09-30 values, and adjusts per-million numbers
# accordingly.

# (Also fills in NaN data, so "cases" is a copy of "corrected_covid_positive".)

import os,json,sys
from os.path import join

tdir="zoemapdata"

#outdir=tdir
outdir="tempd1"

pop={}
with open(join(tdir,"2020-09-30"),'r') as fp:
  dd=json.load(fp)
  for loc in dd:
    pop[loc]=dd[loc]["population"]

for date in sorted(list(os.listdir(tdir))):
  if 1 or date<"2020-09-30":
    with open(join(tdir,date),'r') as fp:
      dd=json.load(fp)
    for loc in dd:
      d=dd[loc]
      d["population"]=pop[loc]
      d["cases"]=d["corrected_covid_positive"]
      d["cases_string"]="%d"%(int(d["cases"]+0.5))
      p=d["cases"]/pop[loc]
      d["cases_pm"]=p*1e6
      d["cases_pm_string"]="%d"%(int(p*1e6+0.5))
      d["percentage"]=p*100
      d["percentage_string"]="%.1f %%"%(p*100)
    with open(join(outdir,date),'w') as fp:
      json.dump(dd,fp,indent=2)
