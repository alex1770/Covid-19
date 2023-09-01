# Convert GISAID metadata file to day counts:
# <date> <non-variant count> <variant count>
# Suitable of binomial regression

import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-v', '--variant',     type=str,default="BA.2.86", help="Name of variant")
#parser.add_argument('-m', '--mutation',    type=str,default="",   help="Key mutations")
parser.add_argument('-f', '--mindate',     default="2019-01-01",  help="Min sample date of sequence")
parser.add_argument('-t', '--maxdate',     default="9999-12-31",  help="Max sample date of sequence")
parser.add_argument('-g', '--gisaid',      action="store_true",   help="Use GISAID data instead of COG-UK data")
args=parser.parse_args()

from stuff import *
fp=sys.stdin
if args.gisaid:
  keys=["Collection date","Lineage","AA Substitutions"]
  sep="\t"
  keymutations=["NSP2_A31D","Spike_L452W","M_T30A","NSP3_N1708S","Spike_I332V","Spike_V445H","Spike_E484K","N_Q229K","Spike_A264D","Spike_V127F","Spike_S50L"]
else:
  keys=["sample_date","lineage","mutations"]
  sep=","
  keymutations=["orf1ab:A211D","S:L452W","M:T30A","orf1ab:N2526S","S:I332V","S:V445H","S:E484K","N:Q229K","S:A264D","S:V127F","S:S50L"]

d={}
for (date,lineage,mutations) in csvrows_it(fp,keys,sep=sep):
  if len(date)!=10 or date<args.mindate or date>args.maxdate: continue
  if date not in d: d[date]=[0,0]
  if lineage=="": continue
  var=int(lineage.split()[0]==args.variant)
  if not var:
    n=len(keymutations)
    k=len([m for m in keymutations if m in mutations])
    if k>n/2: var=1
  d[date][var]+=1

for date in sorted(list(d)):
  print(date,"%6d %6d"%tuple(d[date]))
