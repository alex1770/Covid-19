#!/usr/bin/pypy

from __future__ import print_function,division
# ^ To enable use of pypy when pypy3 isn't available

from stuff import *
import sys,argparse,platform

reflen=29903# Length of reference genome
targetvariant="BA.2.86"
diffthreshold=1

parser=argparse.ArgumentParser()
parser.add_argument('-m', '--metadata',   default="metadata_sorted_from2023-01-01.tsv",  help="Metadata input file")
parser.add_argument('-f', '--mindate',    default="2023-03-01",                          help="Minimum date selected from metadata input")
parser.add_argument('-d', '--decluster',  action="store_true",                           help="Decluster assumes stdin is in fasta format and includes all IDs selected from metadata")
args=parser.parse_args()

if args.decluster and platform.python_implementation()=="CPython": print("Suggest using PyPy for speed\n")

d={}
for name,loc,date,lin in csvrows(args.metadata,["Virus name","Location","Collection date","Lineage"],sep="\t"):
  if len(date)==10 and date>=args.mindate:
    name=name.replace(' ','_')
    locl=loc.split(" / ")
    if len(locl)<2: continue
    country=locl[0].strip()+" / "+locl[1].strip()
    if len(locl)>=3 and (locl[1]=="USA" or locl[1]=="India"): country+=" / "+locl[2].strip()
    if country not in d: d[country]={}
    if date not in d[country]: d[country][date]=[0,0,[]]
    d[country][date][lin[:len(targetvariant)]==targetvariant]+=1
    d[country][date][2].append((name,loc,lin))

ind={}
for country in d:
  numv=sum(x[1] for x in d[country].values())
  if numv==0: continue
  print(country)
  if args.decluster:
    for date in d[country]:
      for (name,loc,lin) in d[country][date][2]:
        ind[name]=[]
  else:
    country_short=country.split(" / ")[-1].strip().replace(" ","_")
    with open("counts/"+country_short,"w") as fp:
      for date in sorted(list(d[country])):
        print(date,"%6d %6d"%(tuple(d[country][date][:2])),file=fp)

# Bits 3210
#      ACGT
# IUPAC codes: R=AG, Y=CT, K=GT, M=AC, S=CG, W=AT, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT
bit={"A":8, "C":4, "G":2, "T":1, "M":12, "R":10, "W":9, "S":6, "Y":5, "K":3, "V":14, "H":13, "D":11, "B":7, "N":15, "-":15}
def declusternumber(seqs):
  nseq=0
  l=[15]*reflen
  for s in seqs:
    if s=="": nseq+=1# Sequence not known, assumed isolated
    else:
      for i in range(reflen): l[i]&=bit[s[i]]
  sc=[]
  for s in seqs:
    if s!="": sc.append([bit[s[i]] for i in range(reflen) if l[i]==0])
  n=len(sc)
  deg=[0]*n
  m=len([x for x in l if x==0])
  for i in range(n-1):
    for j in range(i+1,n):
      diff=0
      for k in range(m): diff+=(sc[i][k]&sc[j][k]==0)
      similar=(diff<=diffthreshold)
      deg[i]+=similar
      deg[j]+=similar
  for d in deg: nseq+=1/(d+1)
  return nseq

notfound=0
wronglength=0

if args.decluster:
  # Collect relevant sequences
  mode=None
  for x in sys.stdin:
    if x[0]=='>':
      x=x.replace(' ','_')
      f=x.find('|')
      if f>=0: id=x[1:f]
      else: id=x.rstrip()[1:]
      if id in ind: mode=ind[id]
      else: mode=None
    elif mode is not None: mode.append(x.rstrip())

  # Concatenate and check sequences
  for id in ind:
    ind[id]="".join(ind[id])
    if ind[id]=="":
      notfound+=1
      #print("Sequence "+id+" not found",file=sys.stderr)
      #raise RuntimeError("Sequence "+id+" not found")
    elif len(ind[id])!=reflen:
      wronglength+=1
      #raise RuntimeError("Sequence %s has length %d; expecting %d for an aligned sequence"%(id,len(ind[id]),reflen))
  print(notfound,"sequence(s) not found in fasta input",file=sys.stderr)
  print(wronglength,"sequence(s) of the wrong length (not aligned)",file=sys.stderr)
  
  for country in d:
    numv=sum(x[1] for x in d[country].values())
    if numv==0: continue
    for date in d[country]:
      subd={}
      d[country][date][:2]=[0,0]
      for (name,loc,lin) in d[country][date][2]:
        isvar=lin[:len(targetvariant)]==targetvariant
        subd.setdefault((loc,isvar),[]).append(ind[name])
      for (loc,isvar) in subd:
        d[country][date][isvar]+=declusternumber(subd[loc,isvar])
    country_short=country.split(" / ")[-1].strip().replace(" ","_")
    with open("counts/"+country_short,"w") as fp:
      for date in sorted(list(d[country])):
        print(date,"%9.2f %9.2f"%(tuple(d[country][date][:2])),file=fp)
 
