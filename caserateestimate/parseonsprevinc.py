from stuff import *

onsdir="ONSsurvey"

def loadonscsv(name):
  fn=os.path.join(onsdir,name+".csv")
  if not os.path.exists(fn): return None
  ret=[]
  with open(fn,"r") as fp:
    header=next(fp)
    for x in fp:
      l=x.strip().split(',')
      ret.append((Date(l[0]),Date(l[1]),float(l[2]),float(l[3]),float(l[4])))
  return ret
  
def saveonscsv(name,data):
  fn=os.path.join(onsdir,name+".csv")
  with open(fn,"w") as fp:
    print("Date start,Date end,number,low,high",file=fp)
    for row in data:
      print("%s,%s,%.0f,%.0f,%.0f"%(row),file=fp)
  print("Wrote ONS data to",fn)

def readsheet(sheet,keys,denom):
  population=56e6
  cols=[None]*len(keys)
  data=[]
  for row in sheet.iterrows():
    l=list(row[1])
    for (i,key) in enumerate(keys):
      for (j,sheetkey) in enumerate(l):
        if type(sheetkey)==str and key.upper()==sheetkey[:len(key)].upper(): cols[i]=j;break
    if any(col==None for col in cols): continue
    try:
      dates=[Date(x) for x in l[cols[0]].split(" to ")]
      rest=tuple([float(l[col])/denom*population for col in cols[1:]])
      data.append((dates[0],dates[1]+1)+rest)
    except:
      pass
  assert len(data)>0
  return data
  
def onsfilenamedate(fn):
  x=fn.lstrip("covid19infectionsurvey").lstrip("datasets")[:8]
  if int(x[2:4])<=12: x=x[4:]+x[2:4]+x[:2]# The odd date appears in DDMMYYYY format
  return Date(x)

# Find a sheet name similar to targ
# Protects against spurious whitespace and case changes which can happen (e.g. the one published on 2021-12-10)
def sheetname(xl,targ):
  l=xl.sheet_names
  if targ in l: return targ
  t1=targ.strip().upper()
  for x in l:
    if x.strip().upper()==t1: return x
  return None

def getonsprevinc(maxdate="9999-12-31"):
  import os
  maxdate=Date(maxdate)
  
  # ONS Coronavirus Infection Survey data from
  # https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
  
  l=[]
  for x in os.listdir(onsdir):
    if x[-5:]==".xlsx":
      date=onsfilenamedate(x)
      if date<=maxdate: l.append((date,x))
  (latestdate,latestfn)=max(l)
  
  onsprev=loadonscsv(latestdate+".prev")
  onsinc=loadonscsv(latestdate+".inc")
  if onsprev!=None and onsinc!=None: return (onsprev,onsinc)
  
  import pandas as pd
  print("Reading prevalence and incidence from",latestfn)
  fn=os.path.join(onsdir,latestfn)
  xl=pd.ExcelFile(fn)
  
  prev=xl.parse(sheetname(xl,'UK summary - positivity'))
  inc=xl.parse(sheetname(xl,'UK summary - incidence'))
  
  keys=["Time period",
        "Estimated average % of the population",# followed by "that had COVID-19" or "testing positive for COVID-19"
        "95% Lower",# followed by "confidence/credible interval" or other thing
        "95% Upper" # followed by "confidence/credible interval" or other thing
        ]
  onsprev=readsheet(prev,keys,100)
  saveonscsv(latestdate+".prev",onsprev)
  
  keys[1]="Estimated COVID-19 incidence rate per 10,000 people per day"
  onsinc=readsheet(inc,keys,10000)
  saveonscsv(latestdate+".inc",onsinc)

  return (onsprev,onsinc)
