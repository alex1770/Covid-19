import sys,os
from stuff import *

c='Europe / United Kingdom'
mindate='2021-01-01'

#VV=['B.1.1.7','B.1.617','B.1.617.1','B.1.617.2','B.1.617.3']
#VV=['B.1.1.7','B.1.617.2','AY.1','AY.2','AY.3','AY.4','AY.5','AY.6','C.1.2']
#VV=['B.1.617.2','AY.1','AY.2','AY.3','AY.4','AY.5','AY.6','AY.4.2']
#VV=['B.1.617.2','AY.4','AY.4.2','AY.43']
#VV=['B.1.617.2','AY.4','AY.4.2','AY.43','B.1.1.529']
#VV=['B.1','B.1.2','B.1.243','B.1.1.7','B.1.1.519','B.1.427','B.1.429']
#VV=['B.1','B.1.2','B.1.243','B.1.1.7','B.1.1.519','B.1.427','B.1.429']
#VV=['B.1.617.2','AY.*','BA.1','BA.1.1','BA.2','BA.3']
Vv=['Spike_F486V', 'NSP8_N118S', 'BA.2*', 'BA.3', 'BA.4*', 'BA.5*', 'Unassigned']

if len(sys.argv)>1: c=sys.argv[1]
if len(sys.argv)>2: mindate=sys.argv[2]
if len(sys.argv)>3: VV=sys.argv[3].split(',')

try:
  t0=os.path.getmtime('metadata.tsv')
except FileNotFoundError:
  t0=-1e30

try:
  t1=os.path.getmtime('metadata_sorted.tsv')
except FileNotFoundError:
  t1=-1e30

if t0<0 and t1<0: raise FileNotFoundError("Could not find GISAID files metadata.tsv or metadata_sorted.tsv")
if t1>=t0:
  infile='metadata_sorted.tsv';inputsorted=True
else:
  infile='metadata.tsv';inputsorted=False

print('#Using input file',infile)

print("#Country/region:",c)
print("#From:",mindate)

print('#Date              All        ',end='')
for v in VV: print(' %11s'%v,end='')
print('     Others')
d={}
for (date,loc,lineage,mutations) in csvrows(infile,['Collection date','Location','Pango lineage','AA Substitutions'],sep='\t'):
  if loc[:len(c)]!=c: continue
  if len(date)!=10 or date[:2]!="20": continue
  if date<mindate:
    if inputsorted: break
    continue
  if date not in d: d[date]=[0]*(len(VV)+1)
  lineage1=lineage+'.'
  for i,pat in enumerate(VV):
    if lineage==pat: break
    if pat[-1]=='*' and lineage1[:len(pat)-1]==pat[:-1]: break
    if '_' in pat and pat in mutations: break
  else:
    i=len(VV)
  d[date][i]+=1

l=sorted(list(d))
for date in l:
  print(date,end='')
  s=sum(d[date])
  print(" %11d        "%s,end='')
  for n in d[date]:
    print(" %11d"%n,end='')
  print("   ",end='')
  #for n in d[date]:
  #  print("  %6.1f%%"%(n/s*100),end='')
  print()

print('#Date              All        ',end='')
for v in VV: print(' %11s'%v,end='')
print('      Others')
