
import csv,time,calendar
import numpy as np
from subprocess import Popen,PIPE

# (file, columnfrom0)
inp=["zoedatapage/symptomnewcases.deconvolve.csv",# Zoe swab new cases
     "zoedatapage/zoeincidence.deconvolve.csv",1),
     ("confirmed.csv",1)]

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

data=[]
lastday=-1
for (fn,col) in inp:
  with open(fn,'r') as fp:
    r=csv.reader(fp)
    headings=next(r)
    loc=headings[col]
    l=[];startday=None
    for x in r:
      if x[col].strip()!='':
        day=datetoday(x[0])
        if startday==None: startday=day
        if day>lastday: lastday=day
        l.append(float(x[col]))
    data.append([startday,l])

loc='London'#alter
    
# Smooth the confirmed case data
l=[]
n=len(data[2][1])
for i in range(n):
  rad=min(3,i,n-1-i)
  s=sum(data[2][1][i-rad:i+rad+1])/(2*rad+1)
  l.append(s)
data[2][1]=l

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

minday=lastday-60
size=2048
graphfn='combined.London.png'
lw=3
po=Popen("gnuplot",shell=True,stdin=PIPE)
p=po.stdin
write('set terminal pngcairo font "sans,%d" size %d,%d'%(8+size//426,size,size*9//16))
write('set bmargin 6;set lmargin 14;set rmargin 8;set tmargin 5')
write('set output "%s"'%graphfn)
write('set key top center')
title="Three estimates of new cases per million per day in "+loc+" based on processed Zoe data and confirmed case counts"
title+="\\nData sources: https://covid.joinzoe.com/data, https://coronavirus.data.gov.uk as published on or before "+daytodate(lastday)
write('set title "%s"'%title)
write('set xdata time')
write('set format x "%Y-%m-%d"')
write('set timefmt "%Y-%m-%d"')
write('set tics scale 3,0.5')
write('set logscale y')
write('set yrange [100:3000]')
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
for (startday,row) in data:
  for (i,x) in enumerate(row):
    if startday+i>=minday: write(daytodate(startday+i),x)
  write("e")
p.close()
po.wait()
print("Written %s"%graphfn)

