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
  cols=[None]*len(keys)
  data=[]
  for row in sheet.iterrows():
    l=list(row[1])
    for (i,key) in enumerate(keys):
      if key in l: cols[i]=l.index(key)
    if any(col==None for col in cols): continue
    try:
      dates=[Date(x) for x in l[cols[0]].split(" to ")]
      rest=tuple([float(l[col])/denom*population for col in cols[1:]])
      data.append((dates[0],dates[1]+1)+rest)
    except:
      pass
  assert len(data)>0
  return data
  

def getonsprevinc():
  import os
  
  # ONS Coronavirus Infection Survey data from
  # https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
  
  latest=max(x for x in os.listdir(onsdir) if x[-5:]==".xlsx")
  
  onsprev=loadonscsv(latest[:8]+"prev")
  onsinc=loadonscsv(latest[:8]+"inc")
  if onsprev!=None and onsinc!=None: return (onsprev,onsinc)
  
  import pandas as pd
  print("Reading prevalence and incidence from",latest)
  fn=os.path.join(onsdir,latest)
  xl=pd.ExcelFile(fn)
  
  population=56e6
  prev=xl.parse('UK summary - positivity')
  inc=xl.parse('UK summary - incidence')
  
  keys=["Time period",
        "Estimated average % of the population testing positive for COVID-19",
        "95% Lower confidence/credible interval",
        "95% Upper confidence/credible interval"
        ]
  onsprev=readsheet(prev,keys,100)
  saveonscsv(latest[:8]+"prev",onsprev)
  
  keys[1]="Estimated COVID-19 incidence rate per 10,000 people per day"
  onsinc=readsheet(inc,keys,10000)
  saveonscsv(latest[:8]+"inc",onsinc)

  return (onsprev,onsinc)
