# Merge metadata files in GISAID format. List files on command line with '-' for standard input.
# Deduplicate, remove entries with unknown date, filter on sample date>=mindate, remove wastewater samples (Human host only), classify BA.2.86, sort into reverse date order, then output to stdout.
# For duplicate entries, later-specified files take precedence, and all files take precedence over stdin.
# 
# Example usage:
#
# tar xf metadata_tsv_2023_08_31.tar.xz metadata.tsv -O | python3 mergegisaid.py -f - 2023-01-01 gisaid_submitted_2023-08-31.tsv gisaid_submitted_2023-09-01.tsv > metadata_current_sorted.tsv
#
# OR
#
# cat metadata_current_sorted.tsv | python3 mergegisaid.py -f 2023-01-01 - gisaid_submitted_2023-08-31.tsv gisaid_submitted_2023-09-01.tsv > temp && mv temp metadata_current_sorted.tsv

import csv,sys,argparse

parser=argparse.ArgumentParser()
parser.add_argument('-f', '--mindate',  default="2019-01-01",  help="Minimum sample date")
parser.add_argument('incrementfilenames',   nargs='*',         help="Names of increment tsv files")
args=parser.parse_args()

mindate=args.mindate

keys=["Virus name","Accession ID","Collection date","Location","Lineage","AA Substitutions"]
output={}

# Mutations that are highly specific for BA.2.86
keymutations=["NSP2_A31D","Spike_L452W","M_T30A","NSP3_N1708S","Spike_I332V","Spike_V445H","Spike_E484K","N_Q229K","Spike_A264D","Spike_V127F","Spike_S50L"]

def classify(lineage,mutations):
  if len([m for m in keymutations if m in mutations])>len(keymutations)/2: return "BA.2.86"
  return lineage

def add(fp,desc):
  reader=csv.reader(fp,delimiter='\t')
  try:
    headings=next(reader)
  except StopIteration:
    print("Ignoring input:",desc,file=sys.stderr)
    return
  cols=[]
  for k in keys:
    if k=="Lineage" and k not in headings: k="Pango lineage"
    if k not in headings: raise RuntimeError("Could not find heading "+k+" in "+desc)
    cols.append(headings.index(k))
  cd=headings.index("Collection date")
  ai=headings.index("Accession ID")# "Accession ID" seems to be more unique than "Virus Name", which can occur with variations
  li=keys.index("Lineage")
  aa=keys.index("AA Substitutions")
  if "Host" in headings: hn=headings.index("Host")
  else: hn=-1
  for row in reader:
    if len(row[cd])==10 and row[cd]>=mindate:# and (hn==-1 or row[hn]=="Human"):
      row1=[row[col] for col in cols]
      row1[li]=classify(row1[li],row1[aa])
      output[row[ai]]=row1

for fn in args.incrementfilenames:
  if fn=='-': add(sys.stdin,"standard input")
  else:
    with open(fn) as fp:
      add(fp,fn)

cd=keys.index("Collection date")
outputv=list(output.values())
outputv.sort(key=lambda row: row[cd], reverse=True)
writer=csv.writer(sys.stdout,delimiter='\t')
writer.writerow(keys)
for row in outputv: writer.writerow(row)
