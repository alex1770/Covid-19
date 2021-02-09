import csv

for y in [19,20]:
  year=str(y)
  
  # From https://geoportal.statistics.gov.uk/datasets/lsoa-2011-to-clinical-commissioning-group-to-stp-to-cancer-alliances-april-2020-lookup-in-england and similar
  fn='LSOA_(2011)_to_Clinical_Commissioning_Groups_to_Sustainability_and_Transformation_Partnerships_(April_20'+year+')_Lookup_in_England.csv'
  
  # Input headings (and similar for other years)
  # FID,LSOA11CD,LSOA11NM,CCG20CD,CCG20CDH,CCG20NM,STP20CD,STP20NM,CAL20CD,CAL20NM,LAD20CD,LAD20NM
  
  # Output headings (LTLA code, LTLA name, STP code, STP name)
  output=['LAD'+year+'CD','LAD'+year+'NM','STP'+year+'CD','STP'+year+'NM']
  
  # Extract columns and uniquify
  with open(fn,'r') as fp:
    r=csv.reader(fp)
    headings=next(r)
    ii=[headings.index(x) for x in output]
    s=set()
    for x in r:
      s.add(tuple(x[i] for i in ii))
  
  with open('LTLAtoSTP'+year,'w') as fp:
    w=csv.writer(fp)
    w.writerow(output)
    for t in sorted(list(s)):
      w.writerow(t)
  
  with open('LAD'+year+'CD','w') as fp:
    t=set(x[0] for x in s)
    for x in sorted(list(t)):
      print(x,file=fp)

  with open('STP'+year+'NM','w') as fp:
    t=set(x[3] for x in s)
    for x in sorted(list(t)):
      print(x,file=fp)
