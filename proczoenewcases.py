import time,calendar,sys
from os.path import join

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

tdir='zoedatapage'

def processnewcases():
  l=[];ok=0;offset=288
  with open(join(tdir,'zoenewcases'),'r') as fp:
    for x in fp:
      y=x.strip().split()
      l.append((y[0],int(y[1])))
  
  days=[datetoday(x[0]) for x in l]
  assert days==list(range(days[0],days[-1]+1))# check contiguous
  nn=[x[1] for x in l]
  n=len(l)
  period=14
  offset=4
  sameweight=0.1
  
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
  x,resid,rank,sing=np.linalg.lstsq(a,b,rcond=-1)
  
  output=[(daytodate(days[i]-offset), "%9.1f"%nn[i], "%9.1f"%x[period+i-1]) for i in range(n)]
  
  import csv
  from subprocess import Popen,PIPE
  
  csvfn=join(tdir,'zoenewcasesdeconvolve.csv')
  with open(csvfn,'w') as fp:
    writer=csv.writer(fp)
    writer.writerow(['Date','Original Total','Deconvolved'])
    for row in output: writer.writerow(row)
    print("Written",csvfn)
  
  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  for (size,graphfn) in [(2560,'zoenewcasesdeconvolve.png'), (1280,'zoenewcasesdeconvolve.small.png')]:
    po=Popen("gnuplot",shell=True,stdin=PIPE)
    p=po.stdin
    write('set terminal pngcairo font "sans,%d" size %d,%d'%(8+size//426,size,size//2))
    write('set bmargin 6;set lmargin 14;set rmargin %d;set tmargin 5'%(size//256-1))
    write('set output "%s"'%graphfn)
    #write('set key left')
    title="Zoe-estimated new cases per day across the UK from app+swab tests"
    title+="\\nData source: https://covid.joinzoe.com/data"
    write('set title "%s"'%title)
    write('set xdata time')
    write('set format x "%Y-%m-%d"')
    write('set timefmt "%Y-%m-%d"')
    write('set tics scale 3,0.5')
    write('set xtics nomirror')
    write('set xtics "2020-01-06", 86400*7')
    write('set xtics rotate by 45 right offset 0.5,0')
    write('set grid xtics ytics lc rgb "#dddddd" lt 1')
    write('set ylabel "New cases per day"')
    s='plot '
    s+='"-" using 1:2 with linespoints lw 3 title "Cases per day over %d-day period (x-axis date is end of period)", '%period
    s+='"-" using 1:2 with linespoints lw 3 title "Soft-deconvolved cases per day"'
    write(s)
    for row in output: write(row[0],row[1])
    write("e")
    for row in output: write(row[0],row[2])
    write("e")
    p.close()
    po.wait()
    print("Written %s"%graphfn)

if __name__=="__main__":
  processnewcases()
