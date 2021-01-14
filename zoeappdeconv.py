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
name='symptomnewcases'
datafile=join(tdir,'zoesymptomprevalence.csv')
#datafile='zoeselected.csv'

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

# Make recovery kernel: kernel[i] = P(recovery time > i+0.5 days)
def getkernel(shape=2.595, scale=4.48):
  kernel=[]
  for i in range(1000000):
    x=gamma.sf(i+0.5,shape,scale=scale)
    kernel.append(x)
    if x<1e-3: break
  return np.array(kernel)

# Returns map: loc -> list (list of (day-of-reporting, est symptomatic cases on this day)
def loadnewcases():
  days=[]
  with open(datafile,'r') as fp:
    reader=csv.reader(fp)
    locs=next(reader)[1:]
    out={loc:[] for loc in locs}
    for row in reader:
      days.append(datetoday(row[0]))
      for (loc,x) in zip(locs,row[1:]):
        out[loc].append(float(x))
  # Interpolate missing entries
  #for loc in locs:
  #  for i in range(len(days)-1):
  #    d0=days[i];d1=days[i+1]
  #    for day in range(d0+1,d1):
  #      poipoi
  #      for (a,b) in zip(out[loc][i],out[loc][i+1]):
  #        out[loc].append((a*(d1-day)+b*(day-d0))/(d1-d0))
  #  #out[loc].sort()# No
  return days,out

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
  days,l=loadnewcases()
  locs=list(l);locs.sort()
  n=len(days)
  
  output={}
  offset=0
  for loc in locs:
    nn=l[loc]
    kernel=getkernel()
    sameweight=1
    # Get x[i] = new cases between day days[i-1]-offset and days[i]-offset
    x=deconvolve(nn,kernel,sameweight)
    output[loc]=[(daytodate(days[i]-offset), x[i]) for i in range(n)]

  import csv
  from subprocess import Popen,PIPE

  dd={loc: dict(output[loc]) for loc in output}# map: location -> {map: date -> value}
  csvfn=join(tdir,name+'.deconvolve.csv')
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

  for (size,graphfn) in [(2560,name+'.deconvolve.png'), (1280,name+'.deconvolve.small.png')]:
    lw=size//1280-1
    po=Popen("gnuplot",shell=True,stdin=PIPE)
    p=po.stdin
    write('set terminal pngcairo font "sans,%d" size %d,%d'%(8+size//426,size,size//2))
    write('set bmargin 6;set lmargin 14;set rmargin %d;set tmargin 5'%(9))#size//256-1))
    write('set output "%s"'%graphfn)
    write('set key bottom right')
    title="Zoe-estimated symptomatic cases from app in regions/nations, and soft-deconvolved estimate of new cases per million per day"
    title+="\\nEnd of graphs are less reliable.  Data source: https://covid.joinzoe.com/data as published "+daytodate(days[-1])
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

if __name__=="__main__":
  processnewcases()
