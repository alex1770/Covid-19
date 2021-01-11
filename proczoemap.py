from collections import defaultdict
import os,json,csv,time,calendar
from subprocess import Popen,PIPE
from os.path import join

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

keys=["respondent", "population", "corrected_covid_positive"]

locs=["Cambridge", "East Cambridgeshire", "South Cambridgeshire", "Barnet", "Haringey", "Epsom and Ewell", "London", "UK"]

zoetrendfn='zoeselected.csv'

tdir='zoemapdata'
def processdata(tdir):
  dates=os.listdir(tdir)
  dates.sort()
  data={loc:[] for loc in locs}
  pop={}
  shifteddates=[]
  with open(zoetrendfn,'w') as fp:
    writer=csv.writer(fp)
    writer.writerow(['Date']+locs)
    for date in dates:
      tot=defaultdict(float)
      totlon=defaultdict(float)
      with open(join(tdir,date),'r') as fp:
        dd=json.load(fp)
        for d in dd.values():
          for x in keys:
            tot[x]+=d[x]
            if d["region"]=="London": totlon[x]+=d[x]
      #shdate=daytodate(datetoday(date)-1)# Go back a day because Zoe values are reported (and timestamped) the day after they occur
      shdate=date# Change to simple recording by publication date. Adjust for lag later in the pipeline.
      row=[shdate]
      shifteddates.append(shdate)
      for loc in locs:
        if loc=="London": src=totlon
        elif loc=="UK": src=tot
        else: src=dd[loc]
        v=src["corrected_covid_positive"]#/src["population"]*1e3
        row.append("%d"%v)
        data[loc].append(v/src["population"]*1e3)
        if date==dates[-1]: pop[loc]=src['population']
      writer.writerow(row)
  print("Written %s"%zoetrendfn)

  # Smooth the small regions in time
  for loc in locs:
    if pop[loc]<2e6:
      n=len(data[loc])
      newdata=[]
      for i in range(n):
        r=min(3,i,n-1-i)
        newdata.append(sum(data[loc][i-r:i+r+1])/(1+2*r))
      data[loc]=newdata

  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  trendfn='zoeselected.png'
  p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
  write('set terminal pngcairo font "sans,13" size 1920,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%trendfn)
  write('set for [i=9:16] linetype i dashtype (20,7)')
  write('set key left')
  #write('set logscale y')
  title="Zoe-estimated active cases per 1000 people, against publication date"
  write('set title "%s"'%title)
  #write('set xlabel "Days since '+desc+perstring+' reached %g'%thr)
  write('set xdata time')
  write('set format x "%Y-%m-%d"')
  write('set timefmt "%Y-%m-%d"')
  write('set tics scale 3,0.5')
  write('set xtics nomirror')
  write('set xtics "2020-08-31", 86400*7')
  write('set xtics rotate by 45 right offset 0.5,0')
  write('set grid xtics ytics lc rgb "#dddddd" lt 1')
  s='plot '
  for loc in locs:
    if s!='plot ': s+=', '
    s+='"-" using 1:2 with lines lw 2 title "%s"'%loc
  write(s)
  for loc in locs:
    for (d,v) in zip(shifteddates,data[loc]): write(d,v)
    write("e")
  p.close()
  print("Written %s"%trendfn)
  

if __name__=="__main__":
  processdata("zoemapdata")
