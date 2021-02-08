import csv

# From https://geoportal.statistics.gov.uk/datasets/lsoa-2011-to-clinical-commissioning-group-to-stp-to-cancer-alliances-april-2020-lookup-in-england
fn='LSOA__2011__to_Clinical_Commissioning_Group_to_STP_to_Cancer_Alliances__April_2020__Lookup_in_England.csv'

# Input headings
# FID,LSOA11CD,LSOA11NM,CCG20CD,CCG20CDH,CCG20NM,STP20CD,STP20NM,CAL20CD,CAL20NM,LAD20CD,LAD20NM

# Output headings (LTLA code, LTLA name, STP code, STP name)
output=['LAD20CD','LAD20NM','STP20CD','STP20NM']

# Extract columns and uniquify
with open(fn,'r') as fp:
  r=csv.reader(fp)
  headings=next(r)
  ii=[headings.index(x) for x in output]
  s=set()
  for x in r:
    s.add(tuple(x[i] for i in ii))

with open('LTLAtoSTP','w') as fp:
  w=csv.writer(fp)
  w.writerow(output)
  for t in sorted(list(s)):
    w.writerow(t)
