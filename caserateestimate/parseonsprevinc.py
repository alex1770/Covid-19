from stuff import *
from math import sqrt

# ONS Coronavirus Infection Survey data from
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
  
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

def loaddailycsv(name):
  fn=os.path.join(onsdir,name+".csv")
  if not os.path.exists(fn): return None
  ret=[]
  with open(fn,"r") as fp:
    header=next(fp)
    for x in fp:
      l=x.strip().split(',')
      ret.append((Date(l[0]),float(l[1]),float(l[2])))
  return ret
  
def savedailycsv(name,data):
  fn=os.path.join(onsdir,name+".csv")
  with open(fn,"w") as fp:
    print("Date,number,sd",file=fp)
    for row in data:
      print("%s,%.0f,%.0f"%(row),file=fp)
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

def getspreadsheetname(maxdate="9999-12-31"):
  """
    Returns (latestdate,latestfn)
  """
  l=[]
  for x in os.listdir(onsdir):
    if x[-5:]==".xlsx":
      date=onsfilenamedate(x)
      if date<=maxdate: l.append((date,x))
  return max(l)
  

def getonsprevinc(maxdate="9999-12-31"):
  import os
  maxdate=Date(maxdate)
  
  (latestdate,latestfn)=getspreadsheetname(maxdate)
  
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

def ONSconsistencycheck(mindate=Date("2022-01-01"),maxdate=Date("9999-12-31")):
  
  import os
  
  (latestdate,latestfn)=getspreadsheetname(maxdate)
  
  import pandas as pd
  print("Reading ONS spreadsheet",latestfn)
  fn=os.path.join(onsdir,latestfn)
  xl=pd.ExcelFile(fn)
  
  dph=dailyprevhist=xl.parse(sheetname(xl,'1l'))

  latest={}
  earliest={}
  for row in dph.iterrows():
    l=list(row[1])
    if len(l)>=4 and type(l[0])==datetime.datetime and all(type(x)==float for x in l[1:4]):
      d=Date(l[0])
      earliest[d]=l[1:4]
      if d not in latest: latest[d]=l[1:4]

  last=max(latest)
  worst=(0,None)
  s0=s2=0
  for x in earliest:
    if x>=mindate and x<last-28:
      v0=earliest[x]
      v1=latest[x]
      sd=(v0[2]-v0[1])/(2*1.96)
      dev=(v1[0]-v0[0])/sd
      print(x,dev)
      if abs(dev)>abs(worst[0]): worst=(dev,x)
      s0+=1;s2+=dev**2
  print()
  print("Checked ONS daily prevalence from",mindate,"for self-consistency")
  print("Worst deviation = %.2f standard deviations at"%worst[0],worst[1])
  print("Correction factor of standard deviations = %.2f"%(sqrt(s2/s0)),"(i.e., narrow ONS confidence intervals by this factor)")

def getdailyprevalence(maxdate=Date("9999-12-31"),correctionfactor=3.69):
  import os
  from itertools import chain

  population=56e6
  
  (latestdate,latestfn)=getspreadsheetname(maxdate)
  
  dailyprev=loaddailycsv(latestdate+"_dailyprev")
  if dailyprev is not None: return dailyprev
  
  import pandas as pd
  print("Reading ONS spreadsheet",latestfn)
  fn=os.path.join(onsdir,latestfn)
  xl=pd.ExcelFile(fn)
  
  dpl=xl.parse(sheetname(xl,'1b'))# Latest daily prevalence
  dph=xl.parse(sheetname(xl,'1l'))# Historic daily prevalence
  
  data={}
  for row in chain(dpl.iterrows(),dph.iterrows()):
    l=list(row[1])
    if len(l)>=4 and type(l[0])==datetime.datetime and all(type(x)==float for x in l[1:4]):
      d=Date(l[0])
      if d not in data: data[d]=(l[1]/100*population,(l[3]-l[2])/(2*1.96)*correctionfactor/100*population)

  dailyprev=[(d,data[d][0],data[d][1]) for d in sorted(list(data))]

  savedailycsv(latestdate+"_dailyprev",dailyprev)
  
  return dailyprev
