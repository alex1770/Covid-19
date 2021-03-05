
import csv,time,calendar
import numpy as np
from subprocess import Popen,PIPE

inp=["zoedatapage/zoeincidence.deconvolve.csv",# Zoe swab/test-based new cases
     "zoedatapage/symptomnewcases.deconvolve.csv",# Zoe app/symptom-based new cases
     "confirmed.csv"]

# For the moment:
# Treat appnewacses as broadly the correct level
# Treat last date of zoeswabnewcases as correct

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def datetofloatday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d-%H-%M%Z')
  return calendar.timegm(t)/86400

def daytofloatdate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d-%H-%M',t)

pop={
  'England':                   55903504,
  'Wales':                     3138631,
  'South East':                8839864,
  'London':                    8886600,
  'Scotland':                  5438100,
  'East of England':           6486032,
  'Midlands':                  10537679,
  'South West':                5601058,
  'North East and Yorkshire':  8560188,
  'North West':                6992083,
  'Northern Ireland':          1881639,
  'United Kingdom':            66361874,
  'East Midlands':             4761918,
  'West Midlands':             5849119,
  'North East':                2632042,
  'Yorkshire and The Humber':  5428656
}

recipes=[
  ('Midlands', ('West Midlands', 'East Midlands')),
  ('North East and Yorkshire', ('Yorkshire and The Humber', 'North East')),
  ('England', ('South East', 'London', 'East of England', 'Midlands', 'South West', 'North West', 'North East and Yorkshire')),
  ('England', ('South East', 'London', 'East of England', 'Midlands', 'South West', 'North West', 'North East', 'Yorkshire and The Humber')),
  ('United Kingdom', ('England', 'Northern Ireland', 'Scotland','Wales')),
]

data=[]
startday=1e20
endday=-1
fnum=0;out={}
for fn in inp:
  out[fnum]={}
  with open(fn,'r') as fp:
    reader=csv.reader(fp)
    headings=next(reader)
    for row in reader:
      day=datetoday(row[0])
      if day<startday: startday=day
      if day>endday: endday=day
      for (loc,num) in zip(headings[1:],row[1:]):
        if num.strip()=='': continue
        if loc not in out[fnum]: out[fnum][loc]={}
        out[fnum][loc][day]=float(num)
  fnum+=1

# Check days are contiguous
for fn in out:
  for loc in out[fn]:
    l=sorted(list(out[fn][loc]))
    assert l==list(range(l[0],l[-1]+1))
  
# Fill in missing entries as far as possible
for fn in out:
  d=out[fn]
  for (x,Y) in recipes:
    if (x not in d) and all(y in d for y in Y):
      d[x]={}
      for day in d[Y[0]]:
        if all(day in d[y] for y in Y):
          d[x][day]=sum(d[y][day] for y in Y)
          
# Smooth the confirmed case data (must be in position 2)
fn=2
for loc in out[fn]:
  n=len(out[fn][loc])
  m=min(out[fn][loc])
  l=[]
  for i in range(n):
    rad=min(3,i,n-1-i)
    s=sum(out[fn][loc][m+j] for j in range(i-rad,i+rad+1))/(2*rad+1)
    l.append(s)
  for i in range(n):
    out[fn][loc][m+i]=l[i]

# Find locations that are in all of the inputs
locs=set(out[0])
for fn in range(1,fnum): locs.intersection_update(out[fn])

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

for loc in locs:
  minday=endday-60
  size=2048
  graphfn='combined.'+loc.replace(' ','_')+'.png'
  lw=3
  po=Popen("gnuplot",shell=True,stdin=PIPE)
  p=po.stdin
  write('set terminal pngcairo font "sans,%d" size %d,%d'%(8+size//426,size,size*9//16))
  write('set bmargin 6;set lmargin 14;set rmargin 8;set tmargin 5')
  write('set output "%s"'%graphfn)
  write('set key top center')
  title="Three estimates of new cases per million per day in "+loc+" based on processed Zoe data and confirmed case counts. Last few values are less reliable."
  title+="\\nData sources: https://covid.joinzoe.com/data, https://coronavirus.data.gov.uk as published on or before "+daytodate(endday)
  write('set title "%s"'%title)
  write('set xdata time')
  write('set format x "%Y-%m-%d"')
  write('set timefmt "%Y-%m-%d"')
  write('set tics scale 3,0.5')
  write('set logscale y')
  write('set yrange [10:3000]')
  #write('set logscale y 3.1622776601683795')
  #write('set logscale y2')
  write('set xtics nomirror')
  #write('set y2tics')
  #write('set y2label "wer"')
  write('set xtics "2020-01-06", 86400*7')
  write('set xtics rotate by 45 right offset 0.5,0')
  write('set grid xtics ytics lc rgb "#dddddd" lt 1')
  write('set ylabel "New cases per million per day (log scale)"')
  s='plot '
  s+=('"-" using 1:2 with linespoints lw %d title "New cases per million per day based on deconvolved Zoe swab data", ')%lw
  s+=('"-" using 1:2 with linespoints lw %d title "New cases per million per day based on deconvolved Zoe symptom data", ')%lw
  s+=('"-" using 1:2 with linespoints lw %d title "New cases per million per day based on confirmed case counts"')%lw
  write(s)
  for fn in range(fnum):
    days=sorted(list(out[fn][loc]))
    for day in days:
      if day>=minday: write(daytodate(day),out[fn][loc][day]/pop[loc]*1e6)
    write("e")
  p.close()
  po.wait()
  print("Written %s"%graphfn)
