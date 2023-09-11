#!/usr/bin/pypy

from __future__ import print_function,division
# ^ To enable use of pypy when pypy3 isn't available

from stuff import *
import sys,argparse,platform,os

reflen=len(refgenome)
diffthreshold=1

parser=argparse.ArgumentParser()
parser.add_argument('-m', '--metadata',   default="metadata_sorted_from2023-03-01.tsv",  help="Metadata input file")
parser.add_argument('-f', '--mindate',    default="2023-03-01",                          help="Minimum date selected from metadata input")
parser.add_argument('-d', '--decluster',  action="store_true",                           help="Decluster assumes stdin is in fasta format and includes all IDs selected from metadata")
parser.add_argument('-s', '--sorted',     action="store_true",                           help="Use this flag to speed up processing if the input is in reverse date order")
parser.add_argument('-b', '--baseline',   default="",                                    help="Baseline variant(s)")
parser.add_argument('-t', '--target',     default="BA.2.86",                             help="Target variant")
args=parser.parse_args()

if args.decluster and platform.python_implementation()=="CPython": print("Suggest using PyPy for speed\n")

# Input is a GISAID tsv or (from UK) a CLIMB csv
gisaidmode=(args.metadata[-4:]==".tsv")

stats={}
for geneloc in genes.values():
  for p in range(geneloc[0],geneloc[1],3):
    c=refgenome[p:p+3]
    x=codontable.get(c,'?')
    if x not in stats: stats[x]={}
    stats[x][c]=stats[x].get(c,0)+1

# Find most likely triple of bases for each AA
# AA2bases[ new AA name, refgenome triple ] = best guess for new triple
AA2bases={}
w0,w1,w2=0.4,0.7,0.75# Weights that work reasonably well
for x in stats:
  y=stats[x]
  s=sum(stats[x].values())
  for old3 in codontable:
    best=(-1e30,)
    for new3 in y:
      # Mutation from old3 -> new3
      score=y[new3]/s-w0*(old3[0]!=new3[0])-w1*(old3[1]!=new3[1])-w2*(old3[2]!=new3[2])
      if score>best[0]: best=(score,new3)
    AA2bases[x,old3]=best[1]

# Get plausible sequence based on mutation and dropout lists.
def conv_climb_metadata_to_sequence(mutations,ambiguities):
  ml=mutations.split('|')
  am=ambiguities.split('|')
  output=list(refgenome)
  for mut in ml:
    (gene,m)=mut.split(':')
    x=m[0]
    p1=int(m[1:-1])
    y=m[-1]
    if gene=='synSNP':
      p=p1-1
      output[p]=y
    else:
      if gene=='orf1ab':
        if p1<=4401: gene='ORF1a'
        else: gene='ORF1b';p1-=4401
      p0=genes[gene]
      p=p0[0]+(p1-1)*3-1
      output[p:p+3]=list(AA2bases[y,refgenome[p:p+3]])

  for ar in am:
    x=[int(p)-1 for p in ar.split('-')]
    if len(x)==1: output[x[0]]='N'
    else: output[x[0]:x[1]+1]=['N']*(x[1]+1-x[0])
    
  return "".join(output)

if gisaidmode:
  keys=["Virus name","Location","Collection date","Lineage"]
  sep="\t"
  source="GISAID"
else:
  if args.decluster:
    keys=["sequence_name","country","sample_date","lineage","mutations","ambiguities"]
  else:
    keys=["sequence_name","country","sample_date","lineage"]
  sep=","
  source="CLIMB (UK)"

ind={}# Virus name -> sequence (alter name)
d={}
for row in csvrows(args.metadata,keys,sep=sep):
  if len(keys)==4:
    name,loc,date,lin=row
  else:
    name,loc,date,lin,mutations,ambiguities=row
  if len(date)!=10: continue
  if date<args.mindate:
    if args.sorted: break
    continue
  if gisaidmode:
    name=name.replace(' ','_')
    locl=loc.split(" / ")
    if len(locl)<2: continue
    country=locl[0].strip()+" / "+locl[1].strip()
    if len(locl)>=3 and locl[1] in ["USA","India","China"]: country+=" / "+locl[2].strip()
  else:
    country="United Kingdom"
    loc="Europe / United Kingdom / "+name.split('/')[0]
    name="hCoV-19/"+name
    if len(keys)==6: ind[name]=conv_climb_metadata_to_sequence(mutations,ambiguities)
  if country not in d: d[country]={}
  if date not in d[country]: d[country][date]=[0,0,[]]
  v0=int(lin[:len(args.baseline)]==args.baseline)
  v1=int(lin[:len(args.target)]==args.target)
  if v0 or v1:
    d[country][date][v1]+=1# Target variant takes priority if both baseline and target match
    d[country][date][2].append((name,loc,lin))

for country in d:
  num=[sum(x[i] for x in d[country].values()) for i in range(2)]
  if num[0]==0 or num[1]==0: continue
  print(country)
  country_short=country.split(" / ")[-1].strip().replace(" ","_")
  with open(os.path.join("counts",country_short),"w") as fp:
    print("# "+source,file=fp)
    for date in sorted(list(d[country])):
      print(date,"%6d %6d"%(tuple(d[country][date][:2])),file=fp)

# Return the number of "sufficiently different" sequences from the set of sequences given as input.
# Do this by making a graph where one sequences is joined to another if it is near-identical (up to diffthreshold bases different, not counting dropout regions)
# then sum up 1/(degree+1) over all nodes. This would be 1 if it is a single clique, and n if it is n distinct sequences, and otherwise something inbetween.
# Bits 3210
#      ACGT
# IUPAC codes: R=AG, Y=CT, K=GT, M=AC, S=CG, W=AT, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT
bit={"A":8, "C":4, "G":2, "T":1, "M":12, "R":10, "W":9, "S":6, "Y":5, "K":3, "V":14, "H":13, "D":11, "B":7, "N":15, "-":15}
def declusternumber(seqs):
  nseq=0
  # As an optimisation, first reduce to the set of locations which can change amongst the input sequences.
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
  if gisaidmode:
    # Collect relevant sequences
    for country in d:
      for date in d[country]:
        for (name,loc,lin) in d[country][date][2]:
          ind[name]=[]
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
    for id in ind:
      ind[id]="".join(ind[id])
  
  # Check sequences
  for id in ind:
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
    num=[sum(x[i] for x in d[country].values()) for i in range(2)]
    if num[0]==0 or num[1]==0: continue
    for date in sorted(list(d[country])):
      #print("Declustering",country,date)
      subd={}
      d[country][date][:2]=[0,0]
      for (name,loc,lin) in d[country][date][2]:
        isvar=lin[:len(args.target)]==args.target
        subd.setdefault((loc,isvar),[]).append(ind[name])
      for (loc,isvar) in subd:
        d[country][date][isvar]+=declusternumber(subd[loc,isvar])
    country_short=country.split(" / ")[-1].strip().replace(" ","_")
    with open(os.path.join("counts_decluster",country_short),"w") as fp:
      print("# "+source,file=fp)
      for date in sorted(list(d[country])):
        print(date,"%9.2f %9.2f"%(tuple(d[country][date][:2])),file=fp)
 
