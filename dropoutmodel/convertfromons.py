# From https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
# save tables 6a and 6b (countries and regions), then do
# python convertfromons.py <table6aname.csv> <table6bname.csv> > ons_ct.csv

import sys,re

# Input headings: (header line split across multiple lines)
#                 Week starting,Percentage of positive tests- by gene,,,,,,,CT Values- all genes,,,,,,,
#                 ,N protein only,ORF1ab only,S protein only,ORF1ab + N protein ,ORF1ab + S protein,N protein + S protein,ORF1ab + N protein + S protein,Mean,10th Percentile,25th Percentile,50th Percentile,75th Percentile,90th Percentile
#
# Input region title
# Example: "Percentage and CT Values of COVID-19 cases, North East",,,,,,,,,,,,,,,
# Recognise by: '^"Percentage and CT Values of COVID-19 cases, (.*)",'
#
# Input data line
# Example: 14 December 2020,13,7,0,22,2,0,57,25.0,14.7,18.0,26.7,31.9,34.8,,
# Recognise by: '^(\w+ \w+ 202[0-9]),([1-9].*)'

# Output csv headings:
# 0             1       2               3       4       5       6       7       8       9       10
# RegionType	Region	Week started	N only	OR only	S only	OR+N	OR+S	N+S	OR+N+S	Mean
#
# 11               12                   13              14              15
# 10th Percentile  25th Percentile	50th Percentile	75th Percentile	90th Percentile

print("RegionType,Region,Week started,N only,OR only,S only,OR+N,OR+S,N+S,OR+N+S,Mean,10th Percentile,25th Percentile,50th Percentile,75th Percentile,90th Percentile")
for (fn,ty) in [(sys.argv[1],"Country"), (sys.argv[2],"EnglandRegion")]:
  with open(fn,"r") as fp:
    for l in fp:
      m=re.search('^"Percentage and CT Values of COVID-19 cases, (.*)",',l)
      if m!=None: region=m.group(1)
      m=re.search('^(\w+ \w+ 202[0-9]),([1-9].*)',l)
      if m!=None:
        print("%s,%s,%s,%s"%(ty,region,m.group(1),m.group(2).rstrip(',')))
