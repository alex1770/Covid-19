import sys,os,argparse
from stuff import *

#VV=['B.1.1.7','B.1.617','B.1.617.1','B.1.617.2','B.1.617.3']
#VV=['B.1.1.7','B.1.617.2','AY.1','AY.2','AY.3','AY.4','AY.5','AY.6','C.1.2']
#VV=['B.1.617.2','AY.1','AY.2','AY.3','AY.4','AY.5','AY.6','AY.4.2']
#VV=['B.1.617.2','AY.4','AY.4.2','AY.43']
#VV=['B.1.617.2','AY.4','AY.4.2','AY.43','B.1.1.529']
#VV=['B.1','B.1.2','B.1.243','B.1.1.7','B.1.1.519','B.1.427','B.1.429']
#VV=['B.1','B.1.2','B.1.243','B.1.1.7','B.1.1.519','B.1.427','B.1.429']
#VV=['B.1.617.2','AY.*','BA.1','BA.1.1','BA.2','BA.3']
#VV=['Spike_F486V-N_P151S', 'NSP8_N118S', 'BA.2*', 'BA.3', 'BA.4*', 'BA.5*', 'Unassigned']

defaultlineages="Spike_F486V-N_P151S,NSP8_N118S,BA.2*,BA.3,BA.4*,BA.5*,Unassigned"

parser=argparse.ArgumentParser()
parser.add_argument('-b', '--maxbad',      type=float, default=0.05,           help="Maximum proportion of Ns allowed")
parser.add_argument('-l', '--lineages',    default=defaultlineages,            help="Comma-separated list of lineages or mutations; AND together with +/- prefixes")
parser.add_argument('-L', '--location',    default="",                         help="Location prefix; AND together with +/- prefixes")
parser.add_argument('-f', '--mindate',     default="2021-01-01",               help="Min sample date of sequence")
parser.add_argument('-i', '--infile',                                          help="Input tsv file (assumed sorted)")
args=parser.parse_args()

VV=args.lineages.split(',')

if args.infile!=None:
  infile=args.infile
  inputsorted=True
else:
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

# "Compile" expression with '+'s and '-'s into logical components
# E.g., Spike_F486V-N_P151S goes to [(True,"Spike_F486V"),(False,"N_P151S")]
def compile(expr):
  cv=[]
  expr='+'+expr
  while expr!="":
    f0=expr.find('+',1)
    f1=expr.find('-',1)
    n=len(expr)
    if f0==-1: f0=n
    if f1==-1: f1=n
    f=min(f0,f1)
    cv.append((expr[0]=='+',expr[1:f]))
    expr=expr[f:]
  return cv


CVV=[compile(v) for v in VV]
CC=compile(args.location)

VV.append("Others")
numv=len(VV)

print('#Using input file',infile)
print("#Country/region:",args.location)
print("#From:",args.mindate)
print("#Max N-Content:",args.maxbad)

print('#Date         All  ',end='')
wid=[]
for v in VV:
  w=max(len(v)+1,6)
  wid.append(w)
  print(' %*s'%(w,v),end='')
print()
d={}
for (date,loc,lineage,mutations,Ncontent) in csvrows(infile,['Collection date','Location','Pango lineage','AA Substitutions','N-Content'],sep='\t'):
  if len(date)!=10 or date[:2]!="20": continue
  if date<args.mindate:
    if inputsorted: break
    continue
  if Ncontent!="" and float(Ncontent)>args.maxbad: continue
  ok=1
  for wanted,place in CC:
    if (loc[:len(place)]==place)!=wanted: ok=0;break
  if not ok: continue
  if date not in d: d[date]=[0]*numv
  lineage1=lineage+'.'
  for i,cv in enumerate(CVV):
    ok=1
    for (wanted,pat) in cv:
      if (lineage==pat or (pat[-1]=='*' and lineage1[:len(pat)]==pat[:-1]+'.') or ('_' in pat and pat in mutations))!=wanted: ok=0;break
    if ok: break
  else:
    i=numv-1
  d[date][i]+=1

l=sorted(list(d))
for date in l:
  print(date,end='')
  s=sum(d[date])
  print(" %6d  "%s,end='')
  for (w,n) in zip(wid,d[date]):
    print(" %*d"%(w,n),end='')
  print("   ",end='')
  #for n in d[date]:
  #  print("  %6.1f%%"%(n/s*100),end='')
  print()

print('#Date         All  ',end='')
for w,v in zip(wid,VV): print(' %*s'%(w,v),end='')
print()
