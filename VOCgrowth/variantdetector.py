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
  print("Found %d mutations of which %d occur at least %d times in %.3fs"%(len(mutcounts),len(num2name),minmutcount,time.clock()-tim0))
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

def getgrowth(linelist,given=set()):
  daycounts={}
  mutdaycounts=[{} for m in range(nmut)]
  growth={}
  for (day,lin,var,muts) in linelist:
    if given.issubset(muts):
      daycounts[day]=daycounts.get(day,0)+1
      for m in muts:
        mutdaycounts[m][day]=mutdaycounts[m].get(day,0)+1
  for mut in range(nmut):
    m=np.zeros([2,2])
    r=np.zeros(2)
    for day in mutdaycounts[mut]:
      x=day-minday
      vm=mutdaycounts[mut][day]
      v0=daycounts[day]
      if v0>0 and vm>0:
        y=log(vm/v0)
        w=1/(1/v0+1/vm)
        m[0,0]+=w
        m[0,1]+=w*x
        m[1,0]+=w*x
        m[1,1]+=w*x*x
        r[0]+=w*y
        r[1]+=w*x*y
    mi=np.linalg.pinv(m)
    c=np.linalg.solve(m,r)
    cv=[mi[0,0],mi[1,1]]# These should be the variances of c[0],c[1]
    growth[mut]=c[1],sqrt(cv[1])
  return growth


print("GH0",time.clock()-tim0)
growth=getgrowth(linelist)
print("GH1",time.clock()-tim0)
  
sd=6
l=[mut for mut in growth if growth[mut][1]>0]
#l.sort(key=lambda x:-(growth[x][0]-sd*growth[x][1]))
l.sort(key=lambda x:-((growth[x][0]-0.00)/growth[x][1]))
nm=0
for mut in l:
  gr=growth[mut]
  (g,gl,gh)=(gr[0],gr[0]-sd*gr[1],gr[0]+sd*gr[1])
  if gl<0: break
  print("%-20s  %6.3f (%6.3f - %6.3f) %6.2f   %8d"%(num2name[mut],g,gl,gh,gr[0]/gr[1],mutcounts[mut]))
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
    if vm>0 and v0>0: print("         %6.3f"%(log(vm/v0)),end='')
    else: print("              -",end='')
    print(" %5d"%vm,end='')
  print()

