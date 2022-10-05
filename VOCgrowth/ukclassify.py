#!/usr/bin/pypy3
# Similar to classify.py. Should unify if have the time.

import sys,os,argparse
from stuff import *
from variantaliases import aliases
from classify import classify, expandlin, contractlin

infile='cog_metadata_sorted.csv';inputsorted=True
cogdate=datetime.datetime.utcfromtimestamp(os.path.getmtime('cog_metadata.csv.gz')).strftime('%Y-%m-%d')

parser=argparse.ArgumentParser()
parser.add_argument('-f',  '--mindate',     default="2022-01-01",  help="Min sample date of sequence")
parser.add_argument('-t',  '--maxdate',     default="9999-12-31",  help="Max sample date of sequence")
args=parser.parse_args()

mindate=Date(args.mindate)
maxdate=Date(args.maxdate)

# Valid lab locations here are UK, England, Northern_Ireland, Scotland and Wales
# NB lab location isn't necessarily the same as sample location
location="UK"

print("Labs:",location)
print("Date range:",mindate,"-",maxdate)
print("Date of COG metadata file:",cogdate)

fn="cog_metadata_furtherclassified."+cogdate+".csv"
fp=open(fn,'w')
print("sequence_name,sample_date,cog_lineage,new_lineage,full_lineage",file=fp)
for (name,date,p2,lin,mutations) in csvrows(infile,['sequence_name','sample_date','is_pillar_2','lineage','mutations']):
  if location!="UK":
    country=name.split('/')[0]
    if country!=location: continue
  if not (len(date)==10 and date[:2]=="20" and date[4]=="-" and date[7]=="-"): continue
  if date>maxdate: continue
  if date<mindate:
    if inputsorted: break
    continue
  mutations='|'+mutations+'|'

  # If the COG-UK lineage is Unassigned or a prefix of tree-rule lineage, then replace it with tree-rule lineage
  if date>="2022-06-01": newlin=classify(mutations,lin)
  else: newlin=lin

  print(f"{name},{date},{lin},{contractlin(newlin)},{expandlin(newlin)}",file=fp)
fp.close()
print("Written",fn)
