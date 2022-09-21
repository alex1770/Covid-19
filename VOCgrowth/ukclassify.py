#!/usr/bin/pypy3
# Similar to classify.py. Should unify if have the time.

import sys,os,argparse
from stuff import *
from variantaliases import aliases
import classifycog

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

def treeclassify(mutations):
  var=classifycog.treeclassify(mutations)

  # Pro tem: manually classify BA.2.75*, because at time of writing (2022-09-01) the autoclassifier is having a hard time on its own (possibly due to bad training data)
  if ("S:G446S" in mutations)+("S:G257S" in mutations)+("S:K147E" in mutations)+("S:G339H" in mutations)>=2: # BA.2.75*
    if "S:D574V" in mutations:
      if "S:R346T" in mutations: var="BL.1"
      else: var="BA.2.75.1"
    elif "S:R346T" in mutations and "S:F486S" in mutations and "S:D1199N" in mutations: var="BA.2.75.2"
    elif "orf1ab:S1221L" in mutations and "orf1ab:V6107I" in mutations: var="BA.2.75.3"
    elif "S:L452R" in mutations: var="BA.2.75.4"
    elif "S:K356T" in mutations: var="BA.2.75.5"
    else: var="BA.2.75"
  # Manually classify BE.1.1.stuff
  if var=="BE.1.1" and "S:K444T" in mutations: var="BE.1.1.1"
  if var=="BE.1.1.1" and "S:N460K" in mutations: var="BQ.1"
  if var=="BQ.1" and "S:R346T" in mutations: var="BQ.1.1"
  
  # Pro tem: add suffix "+S:R346T" for BA.5.1 as that combination doesn't yet have a designated name
  if var=="BA.5.1" and "S:R346T" in mutations: var+="+S:R346T"

  if var=="Other": var="Unassigned"
  
  return var

ecache={}
def expandlin(lin):
  if lin in ecache: return ecache[lin]
  for (short,long) in aliases:
    s=len(short)
    if lin[:s+1]==short+".": ecache[lin]=long+lin[s:];return ecache[lin]
  ecache[lin]=lin
  return lin

ccache={}
def contractlin(lin):
  if lin in ccache: return ccache[lin]
  lin=expandlin(lin)
  for (short,long) in aliases:
    l=len(long)
    if lin[:l+1]==long+".": ccache[lin]=short+lin[l:];return ccache[lin]
  ccache[lin]=lin
  return lin

#fn="cog_metadata_furtherclassified."+cogdate+".csv"
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
  lin_e=expandlin(lin)

  # If the COG-UK lineage is Unassigned or a prefix of tree-rule lineage, then replace it with tree-rule lineage
  if date>="2022-06-01":
    mylin=treeclassify(mutations)
    mylin_e=expandlin(mylin)
    if lin=="Unassigned" or lin_e==mylin_e[:len(lin_e)]: lin_e=mylin_e

  print(f"{name},{date},{lin},{contractlin(lin_e)},{lin_e}",file=fp)
fp.close()
print("Written",fn)
