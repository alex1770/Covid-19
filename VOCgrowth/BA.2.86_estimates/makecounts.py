#!/usr/bin/pypy

from __future__ import print_function,division
# ^ To enable use of pypy when pypy3 isn't available

from stuff import *
import sys,argparse,platform,os

reflen=len(refgenome)

parser=argparse.ArgumentParser()
parser.add_argument('-m', '--metadata',   default="metadata_sorted_from2023-03-01.tsv",  help="Metadata input file")
parser.add_argument('-f', '--mindate',    default="2023-03-01",                          help="Minimum date selected from metadata input")
parser.add_argument('-d', '--decluster',  action="store_true",                           help="Decluster assumes stdin is in fasta format and includes all IDs selected from metadata")
parser.add_argument('-s', '--sorted',     action="store_true",                           help="Use this flag to speed up processing if the input is in reverse date order")
parser.add_argument('-b', '--baseline',   default="",                                    help="Baseline variant(s)")
parser.add_argument('-t', '--target',     default="BA.2.86",                             help="Target variant")
parser.add_argument('-v', '--verbosity',  type=int, default=1,                           help="Verbosity level (0,1,2,...)")
parser.add_argument('-n', '--dirname',    default="counts",                              help="Prefix of directory name where results are stored")
args=parser.parse_args()

if args.verbosity>=1 and args.decluster and platform.python_implementation()=="CPython": print("Suggest using PyPy for speed\n")

# Input is a GISAID tsv or (from UK) a CLIMB csv
gisaidmode=(args.metadata[-4:]==".tsv")

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
  if args.verbosity>=1: print(country)
  country_short=country.split(" / ")[-1].strip().replace(" ","_")
  os.makedirs(args.dirname,exist_ok=True)
  with open(os.path.join(args.dirname,country_short),"w") as fp:
    print("# "+source,file=fp)
    for date in sorted(list(d[country])):
      print(date,"%6d %6d"%(tuple(d[country][date][:2])),file=fp)

notfound=0
wronglength=0

if args.decluster:
  if gisaidmode:
    # Collect relevant sequences
    # To allow for linebreaks in fasta, construct a list of lines for each virus name (one for each name in the metadata)
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
      if args.verbosity>=2: print("Sequence "+id+" not found",file=sys.stderr)
      #raise RuntimeError("Sequence "+id+" not found")
    elif len(ind[id])!=reflen:
      wronglength+=1
      if args.verbosity>=2: print("Sequence "+id+" is of the wrong length (not aligned)",file=sys.stderr)
      #raise RuntimeError("Sequence %s has length %d; expecting %d for an aligned sequence"%(id,len(ind[id]),reflen))
  if args.verbosity>=1:
    print(notfound,"sequence(s) not found in fasta input",file=sys.stderr)
    print(wronglength,"sequence(s) of the wrong length (not aligned)",file=sys.stderr)
    
  for country in d:
    num=[sum(x[i] for x in d[country].values()) for i in range(2)]
    if num[0]==0 or num[1]==0: continue
    for date in sorted(list(d[country])):
      if args.verbosity>=2: print("Declustering",country,date)
      subd={}
      d[country][date][:2]=[0,0]
      for (name,loc,lin) in d[country][date][2]:
        isvar=lin[:len(args.target)]==args.target
        subd.setdefault((loc,isvar),[]).append(ind[name])
      for (loc,isvar) in subd:
        d[country][date][isvar]+=declusternumber(subd[loc,isvar])
    country_short=country.split(" / ")[-1].strip().replace(" ","_")
    os.makedirs(args.dirname+"_decluster",exist_ok=True)
    with open(os.path.join(args.dirname+"_decluster",country_short),"w") as fp:
      print("# "+source,file=fp)
      for date in sorted(list(d[country])):
        print(date,"%9.2f %9.2f"%(tuple(d[country][date][:2])),file=fp)
 
