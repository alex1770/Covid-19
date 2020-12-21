import time,calendar,sys
from os.path import join
import numpy as np

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

tdir='zoedatapage'
name='zoenewcases'
loc='across the UK'
norm=''

# Returns list of (day-of-reporting, average number of new cases over 14 day period ending at day-of-reporting-4)
def loadnewcases():
  l=[]
  with open(join(tdir,name),'r') as fp:
    for x in fp:
      y=x.strip().split()
      l.append((datetoday(y[0]),float(y[1])))
  days=[x[0] for x in l]
  assert days==list(range(days[0],days[-1]+1))# check contiguous
  return l

# Find x[0,...,n+period-2] such that nn[i] ~= x[i]+x[i+1]+...+x[i+period-1], and x[] doesn't jump too much
def deconvolve(nn,period,sameweight):
  n=len(nn)
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
  x=np.maximum(x,0)
  
  return x

def processnewcases(pubdate):
  l=loadnewcases()
  n=len(l)
  
  days=[x[0] for x in l]
  nn=[x[1] for x in l]

  period=14
  offset=4
  sameweight=0.2
  
  # Hidden variables x[0],...,x[n+period-2], where x[period-1+i] represents the new cases on day days[i]-offset
  x=deconvolve(nn,period,sameweight)

  output=[[(daytodate(days[i]-offset-period//2), "%9.1f"%nn[i]) for i in range(n)],
          [(daytodate(days[i]-offset), "%9.1f"%x[period-1+i]) for i in range(n)]]
  
  import csv
  from subprocess import Popen,PIPE

  d0=dict(output[0])
  d1=dict(output[1])
  s=set(d0);s.update(set(d1))
  csvfn=join(tdir,name+'deconvolve.csv')
  with open(csvfn,'w') as fp:
    writer=csv.writer(fp)
    writer.writerow(['Date','Original Total','Deconvolved'])
    for dt in sorted(list(s)):
      row=[dt,d0.get(dt," "*9),d1.get(dt," "*9)]
      writer.writerow(row)
    print("Written",csvfn)
  
  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  for (size,graphfn) in [(2560,name+'deconvolve.png'), (1280,name+'deconvolve.small.png')]:
    lw=size//640-1
    po=Popen("gnuplot",shell=True,stdin=PIPE)
    p=po.stdin
    write('set terminal pngcairo font "sans,%d" size %d,%d'%(8+size//426,size,size//2))
    write('set bmargin 6;set lmargin 14;set rmargin %d;set tmargin 5'%(size//256-1))
    write('set output "%s"'%graphfn)
    write('set key top center')
    title="Zoe-estimated new cases "+norm+"per day "+loc+" from app+swab tests, and soft-deconvolved version of this"
    title+="\\nData source: https://covid.joinzoe.com/data as published "+pubdate
    write('set title "%s"'%title)
    write('set xdata time')
    write('set format x "%Y-%m-%d"')
    write('set timefmt "%Y-%m-%d"')
    write('set tics scale 3,0.5')
    write('set xtics nomirror')
    write('set xtics "2020-01-06", 86400*7')
    write('set xtics rotate by 45 right offset 0.5,0')
    write('set grid xtics ytics lc rgb "#dddddd" lt 1')
    write('set ylabel "New cases '+norm+'per day"')
    s='plot '
    s+=('"-" using 1:2 with linespoints lw %d title "Cases '+norm+'per day over %d-day period (x-axis date is centre of period)", ')%(lw,period)
    s+=('"-" using 1:2 with linespoints lw %d title "Soft-deconvolved cases '+norm+'per day (end of graph less reliable)"')%lw
    write(s)
    for data in output:
      for row in data: write(row[0],row[1])
      write("e")
    p.close()
    po.wait()
    print("Written %s"%graphfn)

if __name__=="__main__":
  if len(sys.argv)>2: name=sys.argv[2]
  if len(sys.argv)>3: loc=sys.argv[3]
  if len(sys.argv)>4: norm=sys.argv[4]
  processnewcases(sys.argv[1])
