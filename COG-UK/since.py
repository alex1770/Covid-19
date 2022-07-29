# Filter COG metadata by minimum date

import sys,os
from stuff import *

mindate='2022-01-01'
if len(sys.argv)>1: mindate=sys.argv[1]

try:
  t0=os.path.getmtime('cog_metadata.csv')
except FileNotFoundError:
  t0=-1e30

try:
  t1=os.path.getmtime('cog_metadata_sorted.csv')
except FileNotFoundError:
  t1=-1e30

if t0<0 and t1<0: raise FileNotFoundError("Could not find COG-UK files cog_metadata.csv or cog_metadata_sorted.csv")
if t1>=t0:
  infile='cog_metadata_sorted.csv';inputsorted=True
else:
  infile='cog_metadata.csv';inputsorted=False

print("Extracting entries from input file",infile,"since",mindate,file=sys.stderr)

with open(infile) as fp: cr=csv.reader(fp);headings=next(cr)
cw=csv.writer(sys.stdout)
cw.writerow(headings)
for (date,row) in csvrows(infile,['sample_date',None]):
  if len(date)!=10 or date[:2]!="20": continue
  if date<mindate:
    if inputsorted: break
    continue
  cw.writerow(row)
