# Make excess death graph using weekly death stats from https://www.mortality.org/Public/STMF/Outputs/stmf.csv
# Can crosscheck with https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/deathsregistrationsummarytables/2019
# and https://www.euromomo.eu/
# and https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/bulletins/excesswintermortalityinenglandandwales/previousReleases
# and https://www.statista.com/statistics/1131428/excess-deaths-in-england-and-wales/

import os,csv,time,calendar
from subprocess import Popen,PIPE

#countrycode='GBRTENW';countryname='England+Wales'
countrycode='FRATNP';countryname='France'
meanfrom=[2015,2010]
displayfrom=2015
update=False

# YYYY-MM-DD -> day number
def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

fn='stmf.csv'
if update or not os.path.isfile(fn):
  if os.path.exists(fn): os.remove(fn)
  Popen("wget https://www.mortality.org/Public/STMF/Outputs/stmf.csv",shell=True).wait()

cn=-1
wd={}
with open(fn,'r') as fp:
  r=csv.reader(fp)
  for row in r:
    if len(row)>0 and row[0]=='CountryCode':
      cn=row.index('CountryCode')
      yr=row.index('Year')
      wk=row.index('Week')
      sx=row.index('Sex')
      dt=row.index('DTotal')
      continue
    if cn>=0 and row[sx]=='b':
      c=row[cn]
      y=int(row[yr])
      w=int(row[wk])
      d=float(row[dt])
      if y>=min(meanfrom): wd.setdefault(c,[]).append([y,w,d])

n=2
sm=[0.2, 0.5, 1, 0.5, 0.2]# Simple smoothing kernel
for c in wd:
  l=[d for (y,w,d) in wd[c]]
  N=len(l)
  for i in range(N):
    s0=s1=0
    for j in range(-n,n+1):
      k=i+j
      if k>=0 and k<N:
        s0+=sm[j+n]
        s1+=sm[j+n]*l[k]
    wd[c][i][2]=s1/s0

# Yearly deaths for each country
if 0:
  for c in wd:
    tot={};orig=set()
    for (y,w,d) in wd[c]:
      orig.add(y)
      if w<=8: y-=1
      tot[y]=tot.get(y,0)+d
    l=sorted(list(orig))
    for y in l:
      if y>=l[0] and y<2020: print(c,y,tot[y])
    print()
  
tot0={y:[0]*52 for y in meanfrom}
tot1={y:[0]*52 for y in meanfrom}
for (y,w,d) in wd[countrycode]:
  if y<2020:# or w<=8:
    for y0 in meanfrom:
      if y>=y0: tot0[y0][w-1]+=1;tot1[y0][w-1]+=d

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

l=[]
for (y,w,d) in wd[countrycode]:
  if y>=displayfrom:
    a=[daytodate(datetoday("%4d-01-01"%y)+7*(w-1)),d]
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
title+=' average'+'s'*(len(meanfrom)>1)+' for corresponding week of year\\n'
title+='Averaging period excludes 2020. Last date: %s. '%(l[-1][0])
title+='Source: www.mortality.org'
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
