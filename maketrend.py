import csv,sys,getdata
from subprocess import Popen,PIPE
if sys.version_info[0]<3: raise SystemExit("Error: requires Python 3")

# Selection of countries to use for the trend graph
#selectcountries=["US","UK","Italy","Spain","France","Germany","Brazil","Austria","S. Korea","Norway","Sweden","Japan","Australia"]
#selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","Canada","S. Korea","Belgium","Japan","France","Brazil"]
selectcountries=["UK","USA","Italy","Sweden","Germany","Spain","S. Korea","Belgium","France","Brazil","Ecuador","Canada","Ireland","Netherlands","Peru"]

# If perhead is True then count deaths per million population instead of absolute deaths
perhead=True

# Base output filename for trend graph
trendfn0="trend"

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

maxdate="0000-00-00"
      
def getprocessedcountrydata(countries,period=7,perhead=True):
  global maxdate
  processedcases={}
  processeddeaths={}
  for country in countries:
    (dates, confirmed, deaths, recovered, active, newc, newd)=getdata.getcountrydata(country,source=source)
    if len(dates)>0 and dates[-1]>maxdate: maxdate=dates[-1]
    if perhead:
      if country not in pop:
        print("Can't find population for country \"%s\" - skipping"%country,file=sys.stderr)
        continue
      denom=pop[country]/1e6
    else:
      denom=1
    for (inp,output) in [(confirmed,processedcases), (deaths,processeddeaths)]:
      n=len(inp)
      for i in range(n):
        if inp[i]/denom>=thr:
          new=[]
          for j in range(i,n):
            new.append((inp[j]-(inp[j-1] if j>0 else 0))/denom)
          m=len(new)
          for i in range(m-1,-1,-1):
            i0=max(i-period,-1)
            new[i]=sum(new[i0+1:i+1])/(i-i0)
          output[country]=new
          break

  return processedcases,processeddeaths
  
# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

if perhead:
  perstring=' per million'
else:
  perstring=''

cases,deaths=getprocessedcountrydata(selectcountries,period=period,perhead=perhead)

for (stats,desc) in [(cases,'cases'), (deaths,'deaths')]:
  countries=sorted(list(stats))
  trendfn=trendfn0+'_'+desc+'.png'
  
  p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
  write('set terminal pngcairo font "sans,13" size 1920,1280')
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
    for v in stats[country]:
      write(v)
    write("e")
  p.close()
  print("Written trend graph to %s"%trendfn)
  
allcountries=getdata.getcountrylist(source=source)
countries=[x for x in allcountries if x in pop and pop[x]>=minpop]
cases,deaths=getprocessedcountrydata(countries,period=period,perhead=perhead)

for (stats,desc) in [(cases,'cases'), (deaths,'deaths')]:
  countries=list(stats)
  countries.sort(key=lambda x:-stats[x][-1])
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
  write('plot "-" using 2:xtic(1) with boxes')
  for country in countries:
    write('"'+country+'"',stats[country][-1])
  write('quit')
  p.close()
  po.wait()
  print("Written bar chart to %s"%barfn)
