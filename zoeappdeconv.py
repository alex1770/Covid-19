import csv,time,calendar,sys
from math import sqrt
from os.path import join
from scipy.stats import gamma
import numpy as np

# P(t) = estimated symptomatic prevalence at time t
# I(t) = number of new symptomatic cases between time t and t+1
# K(dt) = number of recoveries between time dt-0.5 and dt+0.5 since symptom onset
# P(t+1)  = P(t) + I(t) - (# recovered between t and t+1)
#        ~= P(t) + I(t) - (I(t)*K(0) + I(t-1)*K(1) + ...)
# P(t) = I(t-1) - (# rec between t-1 and t)
#      + I(t-2) - (# rec between t-2 and t-1)
#      + ...
#      = (new cases up to time t) - (recovered up to time t)
#      = I(t-1)+I(t-2)+... - (I(t-1)*P(rec<0.5) + I(t-2)*P(rec<1.5) + I(t-3)*P(rec<2.5) + ...)
#      = I(t-1)*P(rec>0.5) + I(t-2)*P(rec>1.5) + I(t-3)*P(rec>2.5) + ...

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

tdir='zoedatapage'
name='appnewcases'
datafile='zoeregions.csv';place=3
#datafile='zoeselected.csv';place=8

# Make recovery kernel: kernel[i] = P(recovery time > i+0.5 days)
def getkernel(shape=2.595, scale=4.48):
  kernel=[]
  for i in range(1000000):
    x=gamma.sf(i+0.5,shape,scale=scale)
    kernel.append(x)
    if x<1e-3: break
  return np.array(kernel)

# Returns list of (day-of-reporting, est symtomatic cases on this day)
def loadnewcases():
  l=[]
  with open(datafile,'r') as fp:
    r=csv.reader(fp)
    headings=next(r)
    loc=headings[place]
    for x in r:
      l.append((datetoday(x[0]),float(x[place])*1000))# Convert from cases/thousand to cases/million
  days=[x[0] for x in l]
  # Interpolate missing entries
  for i in range(len(days)-1):
    d0=days[i];d1=days[i+1]
    for day in range(d0+1,d1):
      l.append((day,(l[i][1]*(d1-day)+l[i+1][1]*(day-d0))/(d1-d0)))
  l.sort()
  return loc,l

# Given nn[0,...,n-1], find x[-nkern+1,...,n-1] such that
# nn[i] ~= x[i]*kernel[0] + x[i-1]*kernel[1] + ... + x[i-nkern+1]*kernel[nkern-1] for i=0,...n-1
# and x[] doesn't jump too much, then return x[0,...,n-1].
# If nn[i]=number infected at day i then x[i] represents new infections between day i-1 and day i.
def deconvolve(nn,kernel,sameweight):
  # Can't have negative indexes so use x0[i]=x[i-nkern+1] for i=0,...,n+nkern-2
  n=len(nn)
  nkern=len(kernel)
  a=np.zeros((n+n+nkern-1,n+nkern-1))
  b=np.zeros(n+n+nkern-1)
  
  # Make hidden variables predict Zoe-estimated symptomatic cases
  for i in range(n):
    for j in range(i,i+nkern): a[i,j]=kernel[nkern-1+i-j]
    b[i]=nn[i]
  
  # Prior on hidden variables not changing too fast
  for i in range(n+nkern-2):
    a[n+i,i]=-sameweight
    a[n+i,i+1]=sameweight
  
  # Infer least squares best hidden variables
  x0,resid,rank,sing=np.linalg.lstsq(a,b,rcond=-1)
  #x=np.maximum(x,0)
  
  return x0[nkern-1:]

# Not actually useful
def geterr(nn,kernel,sameweight):
  # Can't have negative indexes so use x0[i]=x[i-nkern+1] for i=0,...,n+nkern-2
  n=len(nn)
  err=0
  fp=open('temp','w')
  xnow=deconvolve(nn,kernel,sameweight)
  for n0 in range(1,n+1):
    x=deconvolve(nn[:n0],kernel,sameweight)
    print("%12g %12g"%(xnow[n0-1],x[-1]),file=fp)
    err+=(x[-1]-xnow[n0-1])**2
  fp.close()
  return sqrt(err/n)

def processnewcases():
  loc,l=loadnewcases()
  n=len(l)
  
  days=[x[0] for x in l]
  nn=[x[1] for x in l]
  kernel=getkernel()
  sameweight=1
  
  # Get x[i] = new cases between day days[i-1]-offset and days[i]-offset
  x=deconvolve(nn,kernel,sameweight)

  offset=0
#  output=[#[(daytodate(days[i]-offset), "%9.2f"%nn[i]) for i in range(n)],
#          [(daytodate(days[i]-offset-0.5), "%9.3f"%(x[i])) for i in range(n)]]
  
  output=[#[(daytodate(days[i]-offset), "%9.2f"%nn[i]) for i in range(n)],
          [(daytodate(days[i]-offset), "%9.3f"%(x[i])) for i in range(n)]]
  
  import csv
  from subprocess import Popen,PIPE

  dd=[dict(o) for o in output]
  s=set()
  for d in dd: s.update(set(d))
  csvfn=join(tdir,name+'.deconvolve.csv')
  with open(csvfn,'w') as fp:
    writer=csv.writer(fp)
    #writer.writerow(['Date']+['Col%d'%i for i in range(1,len(dd)+1)])
    writer.writerow(['Date']+[loc])# alter
    for dt in sorted(list(s)):
      row=[dt]+[d.get(dt," "*9) for d in dd]
      writer.writerow(row)
    print("Written",csvfn)
  
  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  for (size,graphfn) in [(2560,name+'.deconvolve.png'), (1280,name+'.deconvolve.small.png')]:
    lw=size//640-1
    po=Popen("gnuplot",shell=True,stdin=PIPE)
    p=po.stdin
    write('set terminal pngcairo font "sans,%d" size %d,%d'%(8+size//426,size,size//2))
    write('set bmargin 6;set lmargin 14;set rmargin %d;set tmargin 5'%(9))#size//256-1))
    write('set output "%s"'%graphfn)
    write('set key top center')
    title="Zoe-estimated symptomatic cases from app in %s, and soft-deconvolved estimate of new cases per day"%loc
    title+="\\nData source: https://covid.joinzoe.com/data as published "+daytodate(days[-1])
    write('set title "%s"'%title)
    write('set xdata time')
    write('set format x "%Y-%m-%d"')
    write('set timefmt "%Y-%m-%d"')
    write('set tics scale 3,0.5')
    write('set xtics nomirror')
    write('set xtics "2020-01-06", 86400*7')
    write('set xtics rotate by 45 right offset 0.5,0')
    write('set grid xtics ytics lc rgb "#dddddd" lt 1')
    #write('set ytics nomirror')
    #write('set ytics 0.1, 10')
    write('set ylabel "Cases per million per day"')
    write('set logscale y 10')
    s='plot '
    s+=('"-" using 1:2 with linespoints lw %d title "Derived estimate of new cases per million per day (end of graph less reliable)"')%lw
    write(s)
    for data in output:
      for row in data: write(row[0],row[1])
      write("e")
    p.close()
    po.wait()
    print("Written %s"%graphfn)

if __name__=="__main__":
  processnewcases()
