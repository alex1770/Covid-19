#!/usr/bin/pypy3
import sys,os
from stuff import *
from math import log

mindate="2022-05-01"
mincount=100

lineages=None;numtoplin=None
if len(sys.argv)>1: mindate=sys.argv[1]
if len(sys.argv)>2:
  if ',' in sys.argv[2]: lineages=sys.argv[2].split(',')
  else: numtoplin=int(sys.argv[2])
if len(sys.argv)>3: mincount=int(sys.argv[3])

infile='metadata.tsv';inputsorted=False
t0=os.path.getmtime(infile)
try:
  t1=os.path.getmtime('metadata_sorted.tsv')
  if t1>=t0: infile='metadata_sorted.tsv';inputsorted=True
except:
  pass

print('#Using input file',infile)
print("#From:",mindate)

allm={}
ml=[]
numl={}
for (date,lineage,mutations) in csvrows(infile,['Collection date','Pango lineage',"AA Substitutions"],sep='\t'):
  if len(date)<10: continue
  if date<mindate:
    if inputsorted: break
    continue
  if lineage=="" or lineage=="Unassigned": continue
  for m in mutations[1:-1].split(','): allm[m]=allm.get(m,0)+1
  ml.append((lineage,mutations))
  numl[lineage]=numl.get(lineage,0)+1
  
print("Found",len(ml),"relevant entries since",mindate)
okm=set(m for m in allm if m!="" and allm[m]>=mincount)
print("Found",len(allm),"mutations, of which",len(okm),"have occurred at least",mincount,"times")
if lineages==None: lineages=list(set(lineage for (lineage,mutations) in ml))
lineages.sort(key=lambda l: -numl[l])
if numtoplin!=None: lineages=lineages[:numtoplin]
mml=[]
for (lineage,mutations) in ml:
  if lineage in lineages: i=lineages.index(lineage)
  else: i=len(lineages)
  mml.append((i,set(mutations[1:-1].split(',')).intersection(okm)))
lineages.append("Others")
print("Made reduced list")

def ent(mlist):
  mlist=list(mlist)
  n=len(mlist)
  nl=len(lineages)
  count=[[0]*nl for j in range(1<<n)]
  for (lineage,mutations) in mml:
    i=0
    for (k,m) in enumerate(mlist):
      if m in mutations: i+=1<<k
    # i=0,...,2^n-1, encoding the predictor subset of the specified mlist
    count[i][lineage]+=1
  e=0
  for i in range(1<<n):
    s=float(sum(count[i]))
    e+=sum(x*log(x/s) for x in count[i] if x>0)
  return e

def pr(mlist):
  mlist=list(mlist)
  n=len(mlist)
  nl=len(lineages)
  count=[[0]*nl for j in range(1<<n)]
  maxcol=20
  for (lineage,mutations) in mml:
    i=0
    for (k,m) in enumerate(mlist):
      if m in mutations: i+=1<<k
    # i=0,...,2^n-1, encoding the predictor subset of the specified mlist
    count[i][lineage]+=1
  print(" "*n,end='')
  for j in range(min(nl,maxcol)): print(" %9s"%lineages[j],end='')
  if nl>maxcol: print(" ...",end='')
  print()
  for i in range(1<<n):
    print(''.join(str(i>>k&1) for k in range(n)),end='')
    for j in range(min(nl,maxcol)):
      print(" %9d"%count[i][j],end='')
    print()
  for m in mlist: print(m)

chosen=[]
e0=ent(chosen)
print("Start logprob =",e0)
pr(chosen);print()
while len(chosen)<12:
  best=(e0,)
  for m in okm:
    e=ent(chosen+[m])
    if e>best[0]: best=(e,m)
    #print(m,e)
  print("Best:",best)
  e0=best[0]
  chosen.append(best[1])
  pr(chosen);print()
  sys.stdout.flush()
