import time,calendar,sys,csv
from os.path import join
import numpy as np

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

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

tdir='zoedatapage'
incsv='zoeincidence.csv'
outname='zoeincidence.deconvolve'

def loadnewcases():
  l=[]
  with open(join(tdir,incsv),'r') as fp:
    reader=csv.reader(fp)
    locs=next(reader)[1:]
    out={loc:[] for loc in locs}
    days=[];rowvals=[]
    for row in reader:
      rowval=[float(x) for x in row[1:]]
      day=datetoday(row[0])
      if len(days)>0 and day>days[-1]+1:
        # Interpolate missing entries
        prevval=rowvals[-1]
        for d in range(days[-1]+1,day):
          a,b=d-days[-1],day-d
          vals=[(b*n0+a*n1)/(a+b) for (n0,n1) in zip(prevval,rowval)]
          rowvals.append(vals)
          days.append(d)
      rowvals.append(rowval)
      days.append(day)
    for rowval in rowvals:
      for (loc,x) in zip(locs,rowval):
        out[loc].append(x)
  assert days==list(range(days[0],days[-1]+1))# check contiguous
  return days,out

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
  
  return x[period-1:]

if __name__=="__main__":
  days,l=loadnewcases()
  n=len(days)
  locs=list(l)
  output={}
  
  for loc in locs:
    nn=l[loc]
    period=14
    offset=4
    sameweight=0.2
    
    # Hidden variables x[0],...,x[n+period-2], where x[period-1+i] represents the new cases on day days[i]-offset
    x=deconvolve(nn,period,sameweight)
    output[loc]=[(daytodate(days[i]), x[i]) for i in range(n)]
    
  from subprocess import Popen,PIPE
  dd={loc: dict(output[loc]) for loc in output}# map: location -> {map: date -> value}
  csvfn=join(tdir,outname+'.csv')
  with open(csvfn,'w') as fp:
    writer=csv.writer(fp)
    writer.writerow(['Date']+locs)
    for day in days:
      dt=daytodate(day)
      row=[dt]+[("%.1f"%(dd[loc][dt]) if dt in dd[loc] else " "*9) for loc in locs]
      writer.writerow(row)
    print("Written",csvfn)
  
  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  name=loc.replace(' ','_')
  for (size,graphfn) in [(2560,outname+'.png'), (1280,outname+'.small.png')]:
    lw=size//640-1
    po=Popen("gnuplot",shell=True,stdin=PIPE)
    p=po.stdin
    write('set terminal pngcairo font "sans,%d" size %d,%d'%(8+size//426,size,size//2))
    write('set bmargin 6;set lmargin 14;set rmargin %d;set tmargin 5'%(size//256-1))
    write('set output "%s"'%graphfn)
    write('set key top center')
    title="Soft deconvolved version of Zoe-estimated new cases per million per day from app+swab tests"
    title+="\\nData source: https://covid.joinzoe.com/data as published on "+daytodate(days[-1])
    write('set title "%s"'%title)
    write('set xdata time')
    write('set format x "%Y-%m-%d"')
    write('set timefmt "%Y-%m-%d"')
    write('set tics scale 3,0.5')
    write('set xtics nomirror')
    write('set xtics "2020-01-06", 86400*7')
    write('set xtics rotate by 45 right offset 0.5,0')
    write('set grid xtics ytics lc rgb "#dddddd" lt 1')
    write('set ylabel "New cases per million per day"')
    s='plot '
    s='plot ';first=True
    for loc in locs:
      if not first: s+=", "
      s+=('"-" using 1:2 with linespoints lw %d title "%s"')%(lw,loc)
      first=False
    write(s)
    for loc in locs:
      for row in output[loc]: write(row[0],row[1]/pop[loc]*1e6)
      write("e")
    p.close()
    po.wait()
    print("Written %s"%graphfn)
