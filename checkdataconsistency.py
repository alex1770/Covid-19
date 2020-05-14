import csv,sys
if sys.version_info[0]<3: print("Requires Python 3");sys.exit(1)

# Convert M/D/YY to YYYY-MM-DD
def isofy(d):
  x=d.split('/');assert len(x)==3
  return "20"+x[2]+"-"+("0"+x[0])[-2:]+"-"+("0"+x[1])[-2:]

# List countries
if 0:
  for source in ["worldometer.csv","ecdc.csv","time_series_covid19_confirmed_global.csv","time_series_covid19_deaths_global.csv"]:
    with open(source) as fp:
      r=csv.reader(fp)
      first=1
      for x in r:
        if first: first=0;continue
        print(x[1])

equivnames={}
with open("countrynames") as fp:
  r=csv.reader(fp)
  for x in r:
    if len(x)==0 or x[0][:1]=='#': continue
    for y in x[1:]: equivnames[y]=x[0]

# 2020-03-01 -> 2020-02-29 etc
def yesterday(x):
  (y,m,d)=map(int,x.split('-'))
  d-=1
  if d==0:
    m-=1
    if m==0: y-=1;m=12
    d=[31,28,31,30,31,30,31,31,30,31,30,31][m-1]+(m==2 and y%4==0)
  return "%04d-%02d-%02d"%(y,m,d)
    
d={}
for (source,cols) in [("worldometer.csv",(2,3)),("ecdc.csv",(4,5))]:
  with open(source) as fp:
    r=csv.reader(fp)
    first=1
    for x in r:
      if first: first=0;continue
      date=x[0]
      if source[:4]=="ecdc": date=yesterday(date)
      country=equivnames.get(x[1],x[1])
      key=(country,date)
      if key not in d: d[key]={}
      try:
        d[key][source]=[int(x[cols[0]]),int(x[cols[1]])]
      except:
        pass

dates=None
for ty in ["confirmed","deaths"]:
  with open("time_series_covid19_"+ty+"_global.csv") as fp:
    r=csv.reader(fp)
    first=1
    for x in r:
      if first:
        if dates==None:
          dates=x[4:]
        else:
          assert dates==x[4:]
        first=0
        continue
      if x[0]=="":
        country=equivnames.get(x[1],x[1])
        for (date,vl) in zip(dates,x[4:]):
          date=isofy(date)
          key=(country,date) 
          if key not in d: d[key]={}
          if "jhu" not in d[key]: d[key]["jhu"]=["?","?"]
          d[key]["jhu"][ty=="deaths"]=int(vl)

for key in d:
  (country,date)=key
  ydate=yesterday(date)
  ykey=(country,ydate)
  for ds in d[key]:
    if ykey in d and ds in d[ykey]: prv=d[ykey][ds]
    else: prv=[0,0]
    for ty in [0,1]:
      d[key][ds].append(d[key][ds][ty]-prv[ty])

W='worldometer.csv'
E='ecdc.csv'
J='jhu'
err={}
errt={}
errl=[]
for key in sorted(list(d)):
  (country,date)=key
  if country not in err: err[country]={W:0, E:0, J:0}
  vl=d[key]
  if len(vl)==3:
    for ty in range(2):# range(4) to incorporate newcases, newdeaths
      l=sorted([vl[W][ty],vl[E][ty],vl[J][ty]])
      if l[1]-l[0]<l[2]-l[1]: v=(l[0]+l[1])/2
      else: v=(l[1]+l[2])/2
      for ds in [W,E,J]:
        e=abs((vl[ds][ty]+10)/(v+10)-1)
        errl.append((e,ds,country,date,["cases","deaths","newcases","newdeaths"][ty],vl[ds][ty],int(v)))
        if e>0.05:
          #print("%20s %s %20s  %d  %6d %6d"%(country,date,ds,ty,vl[ds][ty],int(v)))
          err[country][ds]+=1
          errt[ds]=errt.get(ds,0)+1
    #print(key,"-",vl[W],vl[E],vl[J],err)

countries=sorted(list(set(country for (country,date) in d)))
for country in countries:
  if err[country]!={W:0, E:0, J:0}:
    print("%30s"%country,err[country])

print()
print(errt)

print()
print(" Error              Dataset              Country        Date       Type   Given Correct(maybe)")
errl.sort(reverse=True)
for (e,ds,country,date,ty,given,shouldbe) in errl:
  if e<0.01: break
  print("%6.3f %20s %30s  %s  %9s %7d %7d"%(e,ds,country,date,ty,given,shouldbe))
  
