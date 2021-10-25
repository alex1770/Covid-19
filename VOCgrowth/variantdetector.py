import sys,time,os,pickle
from stuff import *
import numpy as np
from math import log,sqrt

mindate='2021-07-01'
maxdate='9999-12-31'

minday=datetoday(mindate)
minmutcount=50# Ignore mutations that have occurred less than this often
cachedir='cogdatacachedir'
datafile='cog_metadata.csv'

tim0=time.clock()
datamtime=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile)).strftime('%Y-%m-%d-%H-%M-%S')
id=datamtime+'_'+mindate+'_'+maxdate+'_'+str(minmutcount)
fn=os.path.join(cachedir,id)
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    (linelist,num2name,name2num,mutcounts,daycounts,mutdaycounts)=pickle.load(fp)
  nmut=len(num2name)
else:
  print("Reading %s"%datafile)
  mutcounts={}
  with open(datafile,'r') as fp:
    for (dt,lin,var,mut) in csvrows_it(fp,['sample_date','lineage','scorpio_call','mutations']):
      if dt<mindate or dt>maxdate: continue# (can't exit early because dates can be out of order)
      day=datetoday(dt)
      muts=mut.split('|')
      for mut in muts:
        mutcounts[mut]=mutcounts.get(mut,0)+1
  num2name=[]
  name2num={}
  for mut in sorted(mutcounts):
    if mutcounts[mut]>=minmutcount:
      name2num[mut]=len(num2name)
      num2name.append(mut)
  print("Found %d mutations of which %d occur at least %d times, in %.3fs"%(len(mutcounts),len(num2name),minmutcount,time.clock()-tim0))
  nmut=len(num2name)
  linelist=[]
  mutcounts=[mutcounts[mut] for mut in num2name]
  daycounts={}
  mutdaycounts=[{} for m in range(nmut)]
  with open(datafile,'r') as fp:
    for (dt,lin,var,mut) in csvrows_it(fp,['sample_date','lineage','scorpio_call','mutations']):
      if dt<mindate or dt>maxdate: continue# (can't exit early because dates can be out of order)
      day=datetoday(dt)
      daycounts[day]=daycounts.get(day,0)+1
      muts=mut.split('|')
      l=[]
      for mut in muts:
        if mut in name2num: l.append(name2num[mut])
      linelist.append([day,lin,var,l])
      for m in l: mutdaycounts[m][day]=mutdaycounts[m].get(day,0)+1
    linelist.sort(reverse=True)# Sort into reverse date order
  os.makedirs(cachedir,exist_ok=True)
  with open(fn,'wb') as fp:
    pickle.dump([linelist,num2name,name2num,mutcounts,daycounts,mutdaycounts],fp)

def getmutday(linelist,minday1=0,maxday1=1000000,given=set(),lineage=None):# Might need a givennot argument too
  daycounts={}
  mutdaycounts=[{} for m in range(nmut)]
  lincounts={}
  for (day,lin,var,muts) in linelist:
    if day>=minday1 and day<maxday1 and given.issubset(muts) and (lineage==None or lin==lineage):
      daycounts[day]=daycounts.get(day,0)+1
      lincounts[lin]=lincounts.get(lin,0)+1
      for m in muts:
        mutdaycounts[m][day]=mutdaycounts[m].get(day,0)+1
  return daycounts,mutdaycounts,lincounts

def getgrowth(daycounts,mutdaycount):
  # log(varcount/backgroundcount) ~ c0+c1*(day-minday) = growth*(day-crossoverday)
  m=np.zeros([2,2])
  r=np.zeros(2)
  s0=syy=0
  tv0=tv1=0
  for day in mutdaycount:
    x=day-minday
    v1=mutdaycount[day]
    v0=daycounts.get(day,0)-v1
    if v0>0 and v1>0:
      y=log(v1/v0)
      w=1/(1/v0+1/v1)
      m[0,0]+=w
      m[0,1]+=w*x
      m[1,0]+=w*x
      m[1,1]+=w*x*x
      r[0]+=w*y
      r[1]+=w*x*y
      s0+=1
      syy+=w*y*y
      tv0+=v0;tv1+=v1
  if m[0,0]<10: return (0,1000),(0,10),(tv0,tv1)#alter
  mi=np.linalg.pinv(m)
  c=np.linalg.solve(m,r)
  cv=[mi[0,0],mi[1,1]]# These should be the variances of c[0],c[1]
  
  # Want sum( w*(c0+c1*x-y)^2 )
  vr = (c[0]**2*m[0,0] + 2*c[0]*c[1]*m[0,1] - 2*c[0]*r[0] + c[1]**2*m[1,1] - 2*c[1]*r[1] + syy)/s0
  #print(vr)
  # Investigate simple correction for overdispersion. Not sure it's right yet.
  # Try crossing number stats

  return (c[0],sqrt(cv[0])),(c[1],sqrt(cv[1])),(tv0,tv1)
  # This form is nicer to interpret (and minday-independent), but will become singular if c[1]=0:
  # return (minday-c[0]/c[1],sqrt(cv[0])/c[1]),(c[1],sqrt(cv[1]))

print("GH0",time.clock()-tim0)
#daycounts,mutdaycounts,lincounts=getmutday(linelist)
daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-08-01'))
#daycounts,mutdaycounts,lincounts=getmutday(linelist,given={name2num['S:Y145H']})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-07-01'),maxday1=datetoday('2021-08-01'))
#daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-07-01'),maxday1=datetoday('2021-08-01'),given={name2num['S:G142D']})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-01-01'),maxday1=datetoday('2021-08-01'),lineage='AY.4')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-01-01'),maxday1=datetoday('2021-08-15'),given={name2num['S:T95I']},lineage='AY.4')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-07-01'),maxday1=datetoday('2021-08-15'),given={name2num['S:T95I']})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-01-01'),maxday1=datetoday('2021-08-01'),lineage='AY.4')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=datetoday('2021-09-01'),maxday1=datetoday('2021-10-07'),given={name2num['S:Y145H']})
print("GH1",time.clock()-tim0)

if 1:
  growth={};tv={}
  for mut in range(nmut):
    gr=getgrowth(daycounts,mutdaycounts[mut])
    growth[mut]=gr[1]
    tv[mut]=gr[2]
  
  sd=4
  l=[mut for mut in growth if growth[mut][0]>0]
  #l.sort(key=lambda x:-(growth[x][0]-sd*growth[x][1]))
  l.sort(key=lambda x:-((growth[x][0]-0.00)/growth[x][1]))
  nm=0
  for mut in l:
    gr=growth[mut]
    (g,gl,gh)=(gr[0],gr[0]-sd*gr[1],gr[0]+sd*gr[1])
    if gl<0: break
    print("%-20s  %6.3f (%6.3f - %6.3f) %6.2f   %7d %7d"%(num2name[mut],g*100,gl*100,gh*100,gr[0]/gr[1],tv[mut][0],tv[mut][1]))
    nm+=1
  
  nmd=min(nm,10)
  print()
  days=set()
  for mut in l[:nmd]: days.update(mutdaycounts[mut])
  print("                All",end='')
  for mut in l[:nmd]: print(" %20s"%num2name[mut],end='')
  print()
  days=sorted(days)
  for day in days:
    v0=daycounts[day]
    print(daytodate(day),"  %6d"%v0,end='')
    for mut in l[:nmd]:
      vm=mutdaycounts[mut].get(day,0)
      v1=v0-vm
      if vm>0 and v1>0: print("         %6.3f"%(log(vm/v1)),end='')
      else: print("              -",end='')
      print(" %5d"%vm,end='')
    print()
  sys.exit(0)

def pr(daycounts,mutdaycounts,muts):
  print()
  days=set()
  for mut in muts: days.update(mutdaycounts[mut])
  print("                All",end='')
  for mut in muts: print(" %20s"%num2name[mut],end='')
  print()
  days=sorted(days)
  for day in days:
    v0=daycounts[day]
    print(daytodate(day),"  %6d"%v0,end='')
    for mut in muts:
      vm=mutdaycounts[mut].get(day,0)
      v1=v0-vm
      if vm>0 and v1>0: print("         %6.3f"%(log(vm/v1)),end='')
      else: print("              -",end='')
      print(" %5d"%vm,end='')
    print()

window=60
step=30
#minday1=minday
minday1=datetoday('2021-09-01')
prevmuts=set()
currentmuts=set()
sd=10
prevdaycounts,prevmutdaycounts,prevlincounts=getmutday(linelist,minday1=minday1,maxday1=minday1+window)

while 1:
  daycounts,mutdaycounts,lincounts=getmutday(linelist,minday1=minday1,maxday1=minday1+window+1000000,given=currentmuts)#alter
  #l=list(lincounts);l.sort(key=lambda lin:-lincounts[lin])
  #print(l[0],sum(lincounts.values()),lincounts[l[0]])
  if daycounts=={}: break
  growth={}
  for mut in range(nmut):
    if mut not in currentmuts:
      growth[mut]=getgrowth(prevdaycounts,mutdaycounts[mut])[1]

  l=[mut for mut in growth if growth[mut][0]>0]
  l.sort(key=lambda x:-((growth[x][0]-0.00)/growth[x][1]))
  #l.sort(key=lambda x:-(growth[x][0]-sd*growth[x][1]))

  if len(l)>0:
    mut=l[0]
    if growth[mut][0]-sd*growth[mut][1]>0:
      print("%s - %s: adding mutation %s with relative growth %.3f (%.3f - %.3f)"%(daytodate(minday1),daytodate(minday1+window),num2name[mut],growth[mut][0],growth[mut][0]-sd*growth[mut][1],growth[mut][0]+sd*growth[mut][1]))
      currentmuts.add(mut)
      continue
  
  print("No new mutations found at %s - %s"%(daytodate(minday1),daytodate(minday1+window)))
  prevdaycounts=daycounts
  minday1+=step
