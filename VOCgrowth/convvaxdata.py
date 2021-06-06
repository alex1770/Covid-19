import sys,csv

# Convert loose format vaccination data csv from LTLA tab of weekly vaccination data from https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/ to a simple standard csv

# "Under xx", "xx-yy", "xx+"
def parseagerange(x):
  if x[:5]=="Under": return "0-%d"%int(x[6:])
  if '-' in x: y=x.split('-');return "%d-%d"%(int(y[0]),int(y[1])+1)
  if x[-1:]=='+': return "%d-150"%int(x[:-1])
  return None

def flattennumber(x): return int(x.replace(',',''))

reader=csv.reader(sys.stdin)
writer=csv.writer(sys.stdout)
start=False
outrows=[]
for row in reader:
  if 'LTLA Code' in row: headings=row
  if any('0-' in x for x in row):
    for (i,x) in enumerate(row):
      while i>=len(headings): headings.append('')
      if row[i]!='': headings[i]=row[i]
    cols=[]
    outputheadings=[]
    for (i,x) in enumerate(headings):
      if x=='LTLA Code': outputheadings.append(x);lcol=i;continue
      a=parseagerange(x)
      if a!=None: outputheadings.append(a);cols.append(i)
    writer.writerow(outputheadings)
    start=True
    continue
  if start and row[lcol][:1]=='E':
    outrows.append([row[lcol]]+[flattennumber(row[i]) for i in cols])
outrows.sort()
for row in outrows:
  writer.writerow(row)
