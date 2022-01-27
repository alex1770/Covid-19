import sys,os,pickle
from stuff import *

datafile='metadata.tsv'
cachedir='gisaidcachedir'

mindate=Date('2021-11-01')
VV=['BA.1','BA.1.1', 'BA.2']

if len(sys.argv)>1: mindate=Date(sys.argv[1])
if len(sys.argv)>2: VV=sys.argv[2].split(',')

numv=len(VV)
datamtime=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile)).strftime('%Y-%m-%d-%H-%M-%S')
print('Datafile: %s at %s'%(datafile,datamtime))
print('Variants:',' '.join(VV))
print()
id='%s_%s_%s'%(datamtime,mindate,'_'.join(VV))
fn=os.path.join(cachedir,id)
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    (data,continents,countries)=pickle.load(fp)
else:
  continents=set()
  countries=set()
  data={}
  print("Assembling variant counts from GISAID file")
  for (date,location,lineage) in csvrows('metadata.tsv',['Collection date','Location','Pango lineage'],sep='\t'):
    if lineage not in VV or len(date)!=10 or date<mindate: continue
    i=VV.index(lineage)
    l=location.split('/')
    continent,country=l[0].strip(),l[1].strip()
    continents.add(continent)
    countries.add(country)
    for loc in ['World',continent,country]:
      if loc not in data: data[loc]={}
      if date not in data[loc]: data[loc][date]=[0]*numv
      data[loc][date][i]+=1
  
  os.makedirs(cachedir,exist_ok=True)
  with open(fn,'wb') as fp:
    pickle.dump((data,continents,countries),fp)

l=list(data)
l.sort(key=lambda x:sum(y[-1] for y in data[x].values()),reverse=True)

n=10
date0=Date(min(min(data[x]) for x in l[:n]))
date1=Date(max(max(data[x]) for x in l[:n]))
header="      Date"
width=2+6*numv
for loc in l[:n]: header+='%*s'%(width,loc[:width-3])
print(header)
for date in Daterange(date0,date1+1):
  print(date,end='')
  for loc in l[:n]:
    print(' |',end='')
    for x in data[loc].get(str(date),[0]*numv):
      print(' %5d'%x,end='')
  print()
print(header)
