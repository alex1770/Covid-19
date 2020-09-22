#/usr/bin/env python3

# Converts from json with strings only to json with strings, floats, ints

import os,ast,json,sys

floatkeys=["long", "lat", "st_areasha", "st_lengths", "corrected_covid_positive", "cases", "cases_pm", "percentage", "discrete_percentage"]
intkeys=["cartodb_id", "objectid", "bng_e", "bng_n", "respondent", "predicted_covid_positive_count", "population", "discrete_cases_pm"]

print("DISABLED");sys.exit(1)

tdir='zoemapdata'
for x in sorted(os.listdir(tdir)):
  fn=os.path.join(tdir,x)
  with open(fn,'r') as fp:
    d=json.load(fp)
  for k0 in d:
    e=d[k0]
    for k in list(e):
      if k in floatkeys:
        try:
          e[k]=float(e[k])
        except:
          print("Removing",x,k0,k)
          del e[k]
      elif k in intkeys:
        try:
          e[k]=int(e[k])
        except:
          print("Removing",x,k0,k)
          del e[k]

  with open(fn,'w') as fp:
    json.dump(d,fp,indent=2)
