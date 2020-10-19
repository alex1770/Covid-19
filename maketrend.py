import csv,sys,getdata
from subprocess import Popen,PIPE
if sys.version_info[0]<3: raise SystemExit("Error: requires Python 3")

# Selection of countries to use for the trend graph
#selectcountries=["US","UK","Italy","Spain","France","Germany","Brazil","Austria","S. Korea","Norway","Sweden","Japan","Australia"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","Canada","S. Korea","Belgium","Japan","France","Brazil"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Brazil","Mexico","Peru","Chile","Armenia"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Brazil","Mexico","Peru","Chile","Armenia","Panama"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Brazil","Mexico","Peru","Chile","Armenia","Panama","Kyrgyzstan","South Africa","Israel"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Brazil","Peru","Armenia","Panama","South Africa","Israel","Colombia"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Brazil","Peru","Argentina","Panama","South Africa","Israel","Colombia"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Brazil","Argentina","Panama","South Africa","Israel","Costa Rica"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Argentina","Israel","Costa Rica","Switzerland","Netherlands","Czech Republic"]
selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Argentina","Israel","Switzerland","Netherlands","Czech Republic","Armenia"]

# If perhead is True then count deaths per million population instead of absolute deaths
perhead=True

# Base output filename for trend graph with start point threshold for each country
trendfn0="trendthr"

# Base output filename for trend graph with start point at startdate
trendfn1="trendsimple"

# Base output filename for bar chart
barfn0="recent"

# Average over last 'period' days
period=7

# Start the "clock" at this many deaths (or deaths per million in perhead=True mode)
thr=1

# Minimum population for bar chart
minpop=1000000

# Max number of countries to show in bar chart
numbar=20

source="worldometer"
#source="ecdc"

equivnames={}
with open("countrynames") as fp:
  r=csv.reader(fp)
  for x in r:
    if len(x)==0 or x[0][:1]=='#': continue
    for y in x[1:]: equivnames[y]=x[0]

selectcountries=[equivnames.get(x,x) for x in selectcountries]

if perhead:
  pop={}
  with open('population') as fp:
    r=csv.reader(fp)
    for x in r:
      if len(x)==0 or x[0][:1]=='#': continue
      pop[equivnames.get(x[0],x[0])]=int(x[1])

# Start date for plot
startdate="2020-02-15"

# YYYY-MM-DD -> day number
def datetoday(s):
  mm=[0,31,59,90,120,151,181,212,243,273,304,334]
  y=int(s[:4])
  m=int(s[5:7])
  d=int(s[8:10])
  return (y-1970)*365+(y-1969)//4+mm[m-1]+(m>=3 and (y&3)==0)+d-1

def daytodate(n):
  mm=[31,28,31,30,31,30,31,31,30,31,30,31]
  y=1970+(n//1461)*4;n%=1461
  if n>=365: n-=365;y+=1
  if n>=365: n-=365;y+=1
  if n>=366: n-=366;y+=1
  m=0
  while 1:
    o=mm[m]+(m==1 and (y&3)==0)
    if n<o: break
    m+=1;n-=o
  return "%4d-%02d-%02d"%(y,m+1,n+1)

startday=datetoday(startdate)

# Process raw country data by normalising, smoothing over period, restricting to desired date range, determining threshold point, and flagging invalid entries
def processdata(countries,data,period=7,perhead=True):
  maxdate="0000-00-00"
  processedcases={}
  processeddeaths={}
  for country in countries:
    dd=data[country]
    if len(dd)==0: continue
    if perhead:
      if country not in pop:
        print("Can't find population for country \"%s\" - skipping"%country,file=sys.stderr)
        continue
      denom=pop[country]/1e6
    else:
      denom=1
    d0=datetoday(dd[0][0])
    d1=datetoday(dd[-1][0]);assert len(dd)==d1-d0+1
    if dd[-1][0]>maxdate: maxdate=dd[-1][0]
    for (inp,output) in [(1,processedcases), (2,processeddeaths)]:
      l=[]
      t0=None
      for d in range(startday,d1+1):
        i=d-d0
        if i<=0:
          l.append('-')
        else:
          v=dd[i][inp]
          if v=='?': l.append('-')
          else:
            if v/denom>=thr and t0==None: t0=d-startday
            i0=max(i-period,0)
            if dd[i0][inp]=='?': l.append('-')
            else:
              x=(dd[i][inp]-(dd[i0][inp] if i0>=0 else 0))/(i-i0)/denom
              if x>0: l.append(x)# arguably x>=0
              else: l.append('-')
      output[country]=(t0,l)
  return processedcases,processeddeaths,maxdate
  
# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

if perhead:
  perstring=' per million'
else:
  perstring=''

data=getdata.getallsimpledata(source=source)
cases,deaths,maxdate=processdata(selectcountries,data,period=period,perhead=perhead)

for (stats,desc) in [(cases,'cases'), (deaths,'deaths')]:
  countries=sorted(list(stats))
  
  trendfn=trendfn0+'_'+desc+'.png'
  p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
  write('set terminal pngcairo font "sans,13" size 2560,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%trendfn)
  write('set for [i=9:16] linetype i dashtype (20,7)')
  write('set key left')
  write('set logscale y')
  title=("Average new "+desc+perstring+" over last %d day%s, starting when total "+desc+" to date"+perstring+" reached %g")%(period,"" if period==1 else "s",thr)
  title+="\\nSelected countries. Source: %s, %s"%(source,maxdate)
  write('set title "%s"'%title)
  write('set xlabel "Days since '+desc+perstring+' reached %g'%thr)
  write('set grid ytics lc rgb "#dddddd" lt 1')
  s='plot '
  for country in countries:
    if s!='plot ': s+=', '
    s+='"-" using 1 with lines lw 2 title "%s"'%country
  write(s)
  for country in countries:
    (t0,vv)=stats[country]
    for v in vv[t0:]: write(v)
    write("e")
  p.close()
  print("Written trend graph to %s"%trendfn)

  trendfn=trendfn1+'_'+desc+'.png'
  p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
  write('set terminal pngcairo font "sans,13" size 2560,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%trendfn)
  write('set for [i=9:16] linetype i dashtype (20,7)')
  write('set key left')
  write('set logscale y')
  title=("Average new "+desc+perstring+" over last %d day%s, aligned by date")%(period,"" if period==1 else "s")
  title+="\\nSelected countries. Source: %s, %s"%(source,maxdate)
  write('set title "%s"'%title)
  #write('set xlabel "Days since '+desc+perstring+' reached %g'%thr)
  write('set xdata time')
  write('set format x "%Y-%m-%d"')
  write('set timefmt "%Y-%m-%d"')
  write('set tics scale 2,0.5')
  write('set xtics "2020-01-06", 604800')#%startdate)# Date labels on Mondays
  write('set xtics rotate by 45 right offset 0.5,0')
  write('set grid xtics ytics lc rgb "#dddddd" lt 1')
  s='plot '
  for country in countries:
    if s!='plot ': s+=', '
    s+='"-" using 1:2 with lines lw 2 title "%s"'%country
  write(s)
  for country in countries:
    (t0,vv)=stats[country]
    for (i,v) in enumerate(vv): write(daytodate(startday+i),v)
    write("e")
  p.close()
  print("Written trend graph to %s"%trendfn)


allcountries=getdata.getcountrylist(source=source)
countries=[x for x in allcountries if x in pop and pop[x]>=minpop]
cases,deaths,maxdate=processdata(countries,data,period=period,perhead=perhead)

for (stats,desc) in [(cases,'cases'), (deaths,'deaths')]:
  countries=[x for x in stats if stats[x][1][-1]!='-']
  countries.sort(key=lambda x:-stats[x][1][-1])
  countries=countries[:numbar]
  barfn=barfn0+'_'+desc+'.png'
  
  po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
  write('set terminal pngcairo font "sans,12" size 1920,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%barfn)
  write('set key off')
  #write('set label rotate')
  write('set xtics nomirror')
  write('set y2tics mirror')
  title="Average number of "+desc+perstring+" over the last %d day%s"%(period,"" if period==1 else "s")
  title+="\\nTop %d amongst countries with population >=1m. Source: %s, %s"%(numbar,source,maxdate)
  write('set title "%s"'%title)
  write('set grid ytics lc rgb "#dddddd" lt 1')
  write('set boxwidth 0.8')
  write('set style fill solid')
  write('set xtics rotate by 20 right offset 1.5,0')
  write('set yrange [0:]')
  write('plot "-" using 2:xtic(1) with boxes')
  for country in countries:
    write('"'+country+'"',stats[country][1][-1])
  write('quit')
  p.close()
  po.wait()
  print("Written bar chart to %s"%barfn)
