# Source: https://data.europa.eu/euodp/en/data/dataset/QrtzdXsI5w26vnr54SIzpQ giving demo_r_mwk_05.tsv
# Explanation: https://ec.europa.eu/eurostat/cache/metadata/en/demomwk_esms.htm
# Ignoring week 99, which appears to contain no data in the examples I've looked at.
# The week number (1-53) is the ISO 8601 week: weeks run Mon-Sun and are assigned to the year in which Thursday falls.
# For the moment just ignore week 53 and (slightly wrongly) pretend that week i is aligned at the same point in the year for any year.

import os,csv,time,calendar
from subprocess import Popen,PIPE

# See https://ec.europa.eu/eurostat/cache/metadata/en/demomwk_esms.htm for country codes
countrycode='UK';countryname='UK'
#countrycode='FR';countryname='France'
#meanfrom=[2015,2010]
meanfrom=[2015]
displayfrom=2020
update=False

assert displayfrom>=min(meanfrom)

# YYYY-MM-DD -> day number
def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

fn='demo_r_mwk_05.tsv'
if update or not os.path.isfile(fn):
  if os.path.exists(fn): os.remove(fn)
  Popen("wget https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/demo_r_mwk_05.tsv.gz -O -|gunzip -c > %s"%fn,shell=True).wait()

wd={}# Process every country even though only using one at the moment
minr={};maxr={}
s=set()
with open(fn,'r') as fp:
  r=csv.reader(fp,delimiter='\t')
  first=True
  for row in r:
    if first:
      assert row[0][:16]=='age,sex,unit,geo'
      # row[1:] is a list of strings like '2012W39 '; weeks are 01, 02, ..., 52, [53,] 99.
      # Construct indexes, ignoring week 99.
      dates=[]
      for i in range(len(row)-1,0,-1):
        y,w=map(int,row[i].strip().split('W'))
        if y>=min(meanfrom) and w<53:# Just ignore week 53 for the moment. Intend to treat it properly later.
          dates.append((i,y,w,daytodate(datetoday("%4d-01-07"%y)+7*(w-1))))
      first=False
    else:
      (age,sex,_,geo)=row[0].split(',')
      if sex=='T':
        if geo not in wd: wd[geo]={};minr[geo]=1000000;maxr[geo]=0
        s.add(age)
        if age=='TOTAL':
          if age not in wd[geo]: wd[geo][age]=[-1]*len(dates)
          for (r,(i,y,w,d)) in enumerate(dates):
            x=row[i].strip()
            if x[-1]=='p' or x[-1]=='e': x=x[:-1].strip()
            if x!=':':
              wd[geo][age][r]=int(x)
              if r<minr[geo]: minr[geo]=r
              if r>=maxr[geo]: maxr[geo]=r+1
#print(s)            

# Specialise to single country, reduce to valid date range and check contiguous
dd={}
r0=minr[countrycode]
r1=maxr[countrycode]
for age in wd[countrycode]:
  dd[age]=wd[countrycode][age][r0:r1]
  assert -1 not in dd[age]
dates=dates[r0:r1]
N=r1-r0

displayfrom=max(displayfrom,dates[0][1])
for i in range(len(meanfrom)):
  meanfrom[i]=max(meanfrom[i],dates[0][1])

n=2
sm=[0.2, 0.5, 1, 0.5, 0.2]# Simple smoothing kernel
for age in dd:
  l=[]
  for i in range(N):
    s0=s1=0
    for j in range(-n,n+1):
      k=i+j
      if k>=0 and k<N:
        s0+=sm[j+n]
        s1+=sm[j+n]*dd[age][k]
    l.append(s1/s0)
  dd[age]=l

tot0={y:[0]*52 for y in meanfrom}
tot1={y:[0]*52 for y in meanfrom}
for ((i,y,w,d),x) in zip(dates,dd['TOTAL']):
  if y<2020:
    for y0 in meanfrom:
      if y>=y0: tot0[y0][w-1]+=1;tot1[y0][w-1]+=x

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

l=[]
for ((i,y,w,d),x) in zip(dates,dd['TOTAL']):
  if y>=displayfrom:
    a=[daytodate(datetoday("%4d-01-07"%y)+7*(w-1)),x]
    for y0 in meanfrom: a.append(tot1[y0][w-1]/tot0[y0][w-1])
    l.append(a)
    
fn=countryname+'_excess.png'
po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
write('set terminal pngcairo font "sans,13" size 1920,1280')
write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
write('set output "%s"'%fn)
write('set key left')
#write('set logscale y')
title="Mortality in %s compared with "%countryname
title+=','.join("%d"%(2020-y) for y in meanfrom)
title+='-year average'+'s'*(len(meanfrom)>1)+' for corresponding week of year\\n'
title+='Averaging period excludes 2020. Last date: %s. '%(l[-1][0])
title+='Source: Eurostat demo\\\_r\\\_mwk\\\_05'
write('set title "%s"'%title)
write('set grid xtics lc rgb "#e0e0e0" lt 1')
write('set grid ytics lc rgb "#e0e0e0" lt 1')
write('set xtics nomirror')
write('set y2tics mirror')
write('set xtics rotate by 45 right offset 0.5,0')
write('set xdata time')
write('set format x "%Y-%m"')
write('set timefmt "%Y-%m-%d"')
s='plot "-" using 1:2 w lines title "Deaths per week"'
for y in meanfrom:
  s+=', "-" using 1:2 w lines title "Deaths per week, %d year average"'%(2020-y)
write(s)
for x in l: write(x[0],x[1])
write("e")
for i in range(len(meanfrom)):
  for x in l: write(x[0],x[2+i])
  write("e")
p.close()
po.wait()
print("Written %s"%fn)
