# Instead of proczoemap_reg.py
import json,csv,time,calendar
from subprocess import Popen,PIPE

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

zoetrendfn='zoeregions.csv'

mapfile='zoemapdata/2020-12-01'# Use any map output to get the area codes and positions/countries
csvfile='zoedatapage/prevalence.csv'
startdate='2020-09-06'

def processdata_reg():

  # Get mapping data
  with open(mapfile,'r') as fp:
    mapd=json.load(fp)
  regions=set()
  loc={}# map from E-number to region
  pop={}# map from region to population
  for d in mapd.values():
    region0=d["region"]
    co=d["country"]
    if co=="England" or co=="Wales": region=region0
    else: region=co
    regions.add(region)
    loc[d["lad16cd"]]=region
    pop[region]=pop.get(region,0)+d["population"]

  # Read bulk csv data
  dat={};enddate=''
  with open(csvfile,'r') as fp:
    reader=csv.reader(fp)
    for row in reader:
      if row[0][0] not in 'NEWS' or row[4]<startdate: continue
      date=row[4][:10]
      if date>enddate: enddate=date
      k=loc[row[0]],date
      dat[k]=dat.get(k,0)+float(row[1])

  startday=datetoday(startdate)
  endday=datetoday(enddate)
  regions=sorted(list(regions))

  # Write regional csv data
  with open(zoetrendfn,'w') as fp:
    writer=csv.writer(fp)
    writer.writerow(['Date']+regions)
    for day in range(startday,endday+1):
      date=daytodate(day)
      writer.writerow([date]+["%.4g"%(dat[reg,date]/pop[reg]*1e3) for reg in regions])
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
  title="Zoe-estimated active cases per 1000 people, against publication date"
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
  for reg in regions:
    if s!='plot ': s+=', '
    s+='"-" using 1:2 with lines lw 2 title "%s"'%reg
  write(s)
  for reg in regions:
    for day in range(startday,endday+1):
      date=daytodate(day)
      write(date,dat[reg,date]/pop[reg]*1e3)
    write("e")
  p.close()
  print("Written %s"%trendfn)
  #print(pop)
  
if __name__=="__main__":
  processdata_reg()
