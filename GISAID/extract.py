import sys
from stuff import *

c=None
mindate='2020-01-01'
var=None

if len(sys.argv)>1: c=sys.argv[1]
if len(sys.argv)>2: mindate=sys.argv[2]
if len(sys.argv)>3: var=sys.argv[3]

print("Country/region:",c,file=sys.stderr)
print("mindate:",mindate,file=sys.stderr)
print("Variant:",var,file=sys.stderr)

with open('metadata.tsv') as fp: cr=csv.reader(fp,delimiter='\t');headings=next(cr)
cw=csv.writer(sys.stdout,delimiter='\t')
cw.writerow(headings)
for (date,loc,lineage,row) in csvrows('metadata.tsv',['Collection date','Location','Pango lineage',None],sep='\t'):
  if c!=None and loc[:len(c)]!=c: continue
  if len(date)==7: date+='-XX'
  if date<mindate: continue
  if var!=None and lineage!=var: continue
  cw.writerow(row)
