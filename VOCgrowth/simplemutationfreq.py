#!/usr/bin/pypy3

from __future__ import print_function,division
# ^ To enable use of pypy when pypy3 isn't available, or when pypy is faster

import sys,os,argparse,platform
from stuff import *
from math import log
from classify import contractlin, expandlin

parser=argparse.ArgumentParser()
parser.add_argument('-b', '--maxbad',      type=float, default=0.05, help="Maximum proportion of Ns allowed")
parser.add_argument('-f', '--mindate',     default="2019-01-01",     help="Min sample date of sequence")
parser.add_argument('-g', '--gisaid',      action="store_true",      help="Use GISAID data instead of COG-UK data")
parser.add_argument('-t', '--maxdate',     default="9999-12-31",     help="Max sample date of sequence")
parser.add_argument('-m', '--mincount',    type=int, default=50,     help="Min count of mutation: only consider mutations which have occurred at least this many times")
parser.add_argument('-M', '--maxleaves',   type=int, default=50,     help="Maximum number of leaves in the tree")
parser.add_argument('-l', '--lineages',                              help="Comma-separated list of lineages to classify (takes precedence over --numtop)")
parser.add_argument('-n', '--numtoplin',   type=int, default=5,      help="Classify the most prevalent 'n' lineages (alternative to specifying target lineages with --lineages)")
parser.add_argument('-p', '--printtree',   type=int, default=20,     help="Print decision tree with this many leaves")
parser.add_argument('-s', '--genomesubset', type=int, default=4,     help="0 = Only consider mutations in a particular hand-picked subset of RBD, 1 = Only consider mutations in receptor-binding domain, 2 = Only consider mutations in spike gene, 3 = Consider mutations in any gene, but not (COG-designated) synSNPs, 4 = Consider synSNPs that are non-synonymous in some overlapping and functional ORF (e.g., A28330G), 5 = Consider any mutation, including all synSNPs")
args=parser.parse_args()

if args.gisaid:
  infile='metadata_sorted.tsv';inputsorted=True
else:
  infile='cog_metadata_sorted.csv';inputsorted=True

def genomesubsetdesc(s):
  if s==0: return "in hand-picked subset of RBD specified by spike locations "+str(handpickedsubset)
  elif s==1: return "in the RBD"
  elif s==2: return "in the spike gene"
  elif s==3: return "not (COG-designated) synSNPs"
  elif s==4: return "not (COG-designated) synSNPs, except those that overlap with a (thought to be) functional ORF, e.g., A28330G"
  else: return "any location in genome"

# List of overlapping ORFs for which there is evidence that they encode
# https://www.sciencedirect.com/science/article/pii/S0042682221000532
# https://virological.org/t/sars-cov-2-dont-ignore-non-canonical-genes/740/2
# 
accessorygenes=set(range(21744,21861)).union(range(25457,25580)).union(range(28284,28575))

# Hand-picked RBD subset from https://twitter.com/CorneliusRoemer/status/1576903120608600064, https://cov-spectrum.org/collections/54?region=Europe
# plus 144, 252, 484 that I added
handpickedsubset=[144,252,346,356,444,445,446,450,452,460,484,486,490,493,494]

def extractint(s):
  i=0
  while i<len(s):
    if s[i].isdigit(): break
    i+=1
  if i==len(s): return -1
  j=i
  while j<len(s):
    if not s[j].isdigit(): break
    j+=1
  return int(s[i:j])

def okmut(m,lev):
  if lev>=5: return True
  if m[:6]=="synSNP":# Implies COG-UK
    if lev<=3: return False
    loc=extractint(m)
    return loc in accessorygenes
  if lev>=3: return True
  if args.gisaid and m[:6]=="Spike_": m="S:"+m[6:]
  if m[:2]!="S:": return False
  if lev==2: return True
  loc=extractint(m)
  if lev==1: return loc>=329 and loc<=521# RBD
  return loc in handpickedsubset

def mutationlist(mutations):
  if args.gisaid: return mutations[1:-1].split(',')
  else: return mutations.split('|')


targnames=[]
targmuts=None
#f=open('BA.2.86.examples')
if args.gisaid: keys=["Virus name","AA Substitutions"];sep='\t'
else: keys=["sequence_name","mutations"];sep=','
for (name,mutations) in csvrows_it(sys.stdin,keys,sep=sep):
  targnames.append(name)
  muts=[]
  muts=set(mutationlist(mutations))
  if targmuts is None: targmuts=muts
  else: targmuts.intersection_update(muts)
targnames=set(targnames)

print("Found %d target genome(s)"%len(targnames))


print("Reading sequence metadata")
sys.stdout.flush()
allm={mut:0 for mut in targmuts}
t0=t1=0
if args.gisaid: keys=["Virus name","Collection date","Pango lineage","AA Substitutions","N-Content"];sep='\t'
else: keys=["sequence_name","sample_date","usher_lineage","mutations","ambiguities"];sep=','
for (name,date,lineage,mutations,Ncontent) in csvrows(infile,keys,sep=sep):
  if len(date)<10: continue
  if date>args.maxdate: continue
  if date<args.mindate:
    if inputsorted: break
    continue
  if name in targnames: continue
  if lineage=="" or lineage=="Unassigned": continue
  if args.gisaid:
    if Ncontent!="": bad=float(Ncontent)
    else: bad=0
  else:
    t=0
    for x in Ncontent.split('|'):
      y=x.split('-')
      if len(y)==1: t+=1
      else: t+=int(y[1])-int(y[0])+1
    bad=t/29903
  t0+=1
  if bad>args.maxbad: t1+=1;continue
  for m in mutationlist(mutations): allm[m]=allm.get(m,0)+1

print("Found",t0-t1,"relevant entries since",args.mindate)
print("Discarded",t1,"from",t0,"(%.1f%%) due to bad coverage"%(t1/t0*100))

l=sorted(list(targmuts))
lev={mut:min(s for s in range(1,6) if okmut(mut,s)) for mut in targmuts}
l.sort(key=lambda mut:allm[mut]+0.1*lev[mut])
for mut in l:
  if allm[mut]>=10000: break
  print("%-15s %d  %8d"%(mut,lev[mut],allm[mut]))
print()
bestmuts=set(l[:13])

print("Rereading sequence metadata looking for other examples")
sys.stdout.flush()
allm={mut:0 for mut in targmuts}
if args.gisaid: keys=["Virus name","Collection date","Pango lineage","AA Substitutions","N-Content"];sep='\t'
else: keys=["sequence_name","sample_date","usher_lineage","mutations","ambiguities"];sep=','
scorelist=[]
for (name,date,lineage,mutations,Ncontent) in csvrows(infile,keys,sep=sep):
  if len(date)<10: continue
  if date>args.maxdate: continue
  if date<args.mindate:
    if inputsorted: break
    continue
  #if name in targnames: continue
  if lineage=="" or lineage=="Unassigned": continue
  score=0
  for m in mutationlist(mutations):
    if m in bestmuts: score+=1
  scorelist.append((score,name,date))

scorelist.sort(reverse=True)
for (score,name,date) in scorelist[:10]:
  print(date,"%-50s"%name,"%3d"%score,name in targnames)
  
