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

Gretna=(55.006, -3.059)
Berwick=(55.764, -1.999)
Cairngaan=(54.640,-4.917)
Portnahaven=(55.668,-6.529)
def country(lat,long):
  if (lat-Gretna[0])*(Berwick[1]-Gretna[1])-(long-Gretna[1])*(Berwick[0]-Gretna[0])<0: return "England+Wales"
  if (lat-Cairngaan[0])*(Portnahaven[1]-Cairngaan[1])-(long-Cairngaan[1])*(Portnahaven[0]-Cairngaan[0])>0: return "Northern Ireland"
  return "Scotland"

zoetrendfn='zoeregions.csv'

tdir='zoemapdata'
def processdata_reg(tdir):
  dates=os.listdir(tdir)
  dates.sort()
  shifteddates=[]
  sortedregions=None
  with open(zoetrendfn,'w') as fp:
    writer=csv.writer(fp)
    #writer.writerow(['Date']+locs)
    for date in dates:
      tot=defaultdict(float)
      totlon=defaultdict(float)
      totreg={}
      with open(join(tdir,date),'r') as fp:
        dd=json.load(fp)
        for d in dd.values():
          region0=d["region"]
          co=country(d["lat"],d["long"])
          if co=="England+Wales": region=region0
          else: region=co
          if region not in totreg: totreg[region]=defaultdict(float)
          for x in keys:
            tot[x]+=d[x]
            totreg[region][x]+=d[x]
      shdate=daytodate(datetoday(date)-1)# Go back a day because Zoe values are reported (and timestamped) the day after they occur
      sr=sorted(list(totreg))
      if sortedregions==None:
        sortedregions=sr
        writer.writerow(['Date']+sortedregions)
        data={reg:[] for reg in sortedregions}
      else:
        assert sortedregions==sr
      row=[shdate]
      shifteddates.append(shdate)
      for reg in sortedregions:
        src=totreg[reg]
        v=src["corrected_covid_positive"]/src["population"]*1e3
        row.append("%.4g"%v)
        data[reg].append(v)
      writer.writerow(row)
  print("Written %s"%zoetrendfn)

  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  trendfn='zoeregions.png'
  p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
  write('set terminal pngcairo font "sans,13" size 1920,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%trendfn)
  write('set for [i=9:16] linetype i dashtype (20,7)')
  write('set key left')
  #write('set logscale y')
  title="Zoe-estimated active cases per 1000 people"
  write('set title "%s"'%title)
  #write('set xlabel "Days since '+desc+perstring+' reached %g'%thr)
  write('set xdata time')
  write('set format x "%Y-%m-%d"')
  write('set timefmt "%Y-%m-%d"')
  write('set grid xtics ytics lc rgb "#dddddd" lt 1')
  write('set tics scale 3,0.5')
  write('set xtics nomirror')
  write('set xtics "2020-08-31", 86400*7')
  write('set xtics rotate by 45 right offset 0.5,0')
  s='plot '
  for reg in sortedregions:
    if s!='plot ': s+=', '
    s+='"-" using 1:2 with lines lw 2 title "%s"'%reg
  write(s)
  for reg in sortedregions:
    for (d,v) in zip(shifteddates,data[reg]): write(d,v)
    write("e")
  p.close()
  print("Written %s"%trendfn)
  

if __name__=="__main__":
  processdata_reg("zoemapdata")
