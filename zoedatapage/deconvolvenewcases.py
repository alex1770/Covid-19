import time,calendar

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)


l=[];ok=0;offset=288
with open('zoenewcases','r') as fp:
  for x in fp:
    if x[:5]=='START':
      headings=x.strip()
      graphcol=headings.find('2020-11-15')+10;assert graphcol>=10
      datacol=headings.find('datapage')+8;assert datacol>=8
      ok=1
      continue
    if ok:
      if x[:3]=='END': break
      i=x.rfind(' ',0,graphcol)+1
      #if i>graphcol-10 and x[i].isdigit(): v=int(x[i:graphcol])+offset# Prefer graph value, but fall back to datapage value
      if x[:10]<'2020-09-27': v=int(x[i:graphcol])+offset# Prefer graph values up to 2020-09-27
      else:
        i=x.rfind(' ',0,datacol)+1;assert i>0 and x[i].isdigit()
        v=int(x[i:datacol])
      l.append((x[:10],v))

days=[datetoday(x[0]) for x in l]
assert days==list(range(days[0],days[-1]+1))# check contiguous
nn=[x[1] for x in l]
n=len(l)
period=14
offset=4
sameweight=0.1

if 0:
  # Simple deconvolution
  m=[]
  for i in range(1,n):
    s=0
    for j in range(i,0,-period):
      s+=nn[j]-nn[j-1]
    m.append(s*period)

  # Check
  o=[]
  for i in range(n):
    o.append(sum(m[max(i-period,0):i])/period+nn[0])
    assert o[i]==nn[i]
    
  for i in range(n-1):
    print(daytodate(days[i]+1-offset),m[i],nn[i+1])

import numpy as np
# Hidden variables x[0],...,x[n+period-2], where x[period-1+i] represents the new cases at timestep i (reported on day days[i], correponding to days[i]-offset)
a=np.zeros((n+n+period-1,n+period-1))
b=np.zeros(n+n+period-1)

# Make hidden variables predict Zoe-reported averages
for i in range(n):
  for j in range(i,i+period): a[i,j]=1/period
  b[i]=nn[i]

# Prior on hidden variables not changing too fast
for i in range(n+period-2):
  a[n+i,i]=-sameweight
  a[n+i,i+1]=sameweight

# Infer least squares best hidden variables
x,resid,rank,sing=np.linalg.lstsq(a,b,rcond=None)

output=[(daytodate(days[i]-offset), "%9.1f"%nn[i], "%9.1f"%x[period+i-1]) for i in range(n)]

import csv
from subprocess import Popen,PIPE

csvfn='zoenewcasesdeconvolve.csv'
with open(csvfn,'w') as fp:
  writer=csv.writer(fp)
  writer.writerow(['Date','Original Total','Deconvolved'])
  for row in output: writer.writerow(row)
  print("Written",csvfn)

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

graphfn='zoenewcasesdeconvolve.png'
p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
write('set terminal pngcairo font "sans,13" size 2560,1280')
write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
write('set output "%s"'%graphfn)
#write('set for [i=9:16] linetype i dashtype (20,7)')
write('set key left')
title="Zoe-estimated new cases per day from app+swab tests"
write('set title "%s"'%title)
write('set xdata time')
write('set format x "%Y-%m-%d"')
write('set timefmt "%Y-%m-%d"')
write('set tics scale 3,0.5')
write('set xtics nomirror')
write('set xtics "2020-01-06", 86400*7')
write('set xtics rotate by 45 right offset 0.5,0')
write('set grid xtics ytics lc rgb "#dddddd" lt 1')
s='plot '
s+='"-" using 1:2 with linespoints lw 3 title "Cases per day over %d-day period", '%period
s+='"-" using 1:2 with linespoints lw 3 title "Deconvolved cases per day"'
write(s)
for row in output: write(row[0],row[1])
write("e")
for row in output: write(row[0],row[2])
write("e")
p.close()
print("Written %s"%graphfn)
