#/usr/bin/env python3

# Converts from old limited format to json format

import os,ast,json,sys

print("DISABLED");sys.exit(1)

fn=os.path.join('2020-09-12')
with open(fn,'r') as fp:
  staticdict=json.load(fp)

outd={}
fn=sys.argv[1]
sillydate='05:00:00 '+fn[-2:]+'-'+fn[-5:-3]+'-'+fn[-10:-6]
with open(fn,'r') as fp:
  for x in fp:
    d=ast.literal_eval(x)
    loc=d['lad16nm']
    for x in ["objectid", "lad16cd", "lad16nmw", "bng_e", "bng_n", "long", "lat", "st_areasha", "st_lengths", "population", "country", "region"]:
      d[x]=staticdict[loc][x]
    if "data_status" not in d: d["data_status"]=sillydate
    
    y="cases"
    z=y+"_string"
    if y in d and d[y][0].isdigit() and z not in d: d[z]=str(int(float(d[y])+0.5))
    if z in d and y not in d: d[y]=d[z]

    if "cases" in d and d["cases"][0].isdigit() and "population" in d and "cases_pm" not in d:
      d["cases_pm"]="%.3f"%(float(d["cases"])/float(d["population"])*1e6)
    
    y="cases_pm"
    z=y+"_string"
    if y in d and d[y][0].isdigit() and z not in d: d[z]=str(int(float(d[y])+0.5))
    if z in d and y not in d: d[y]=d[z]

    if "corrected_covid_positive" not in d and "cases_string" in d and d["cases_string"][0].isdigit():
      d["corrected_covid_positive"]=d["cases_string"]
    outd[loc]=d

#json.dump(outd,sys.stdout,indent=2)
with open(fn,'w') as fp:
  json.dump(outd,fp,indent=2)

