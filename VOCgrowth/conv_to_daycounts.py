# Convert GISAID metadata file to day counts:
# <date> <non-variant count> <variant count>
# Suitable of binomial regression

import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-v', '--variant',     type=str,default="BA.2.86", help="Name of variant")
parser.add_argument('-f', '--mindate',     default="2019-01-01",  help="Min sample date of sequence")
parser.add_argument('-t', '--maxdate',     default="9999-12-31",  help="Max sample date of sequence")
args=parser.parse_args()

from stuff import *
fp=sys.stdin#open("gisaid_hcov-19_2023_08_22_15.tsv")
d={}
for (date,lineage) in csvrows_it(fp,["Collection date","Lineage"],sep="\t"):
  if len(date)!=10 or date<args.mindate or date>args.maxdate: continue
  if date not in d: d[date]=[0,0]
  if lineage=="": continue
  var=int(lineage.split()[0]==args.variant)
  d[date][var]+=1

for date in sorted(list(d)):
  print(date,"%6d %6d"%tuple(d[date]))
