from collections import defaultdict
import os,json,csv,time,calendar
import pandas as pd
from matplotlib import pyplot
from os.path import join

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

keys=["respondent", "population", "corrected_covid_positive"]

locs=["Cambridge", "East Cambridgeshire", "South Cambridgeshire", "Barnet", "Haringey", "Epsom and Ewell", "London", "UK"]

selectedfn='zoeselected.csv'

tdir='zoemapdata'
if 1:#def processdata(tdir):
  l=os.listdir(tdir)
  l.sort()
  with open(selectedfn,'w') as fp:
    writer=csv.writer(fp)
    writer.writerow(['Date']+locs)
    for date in l:
      tot=defaultdict(float)
      totlon=defaultdict(float)
      with open(join(tdir,date),'r') as fp:
        dd=json.load(fp)
        for d in dd.values():
          for x in keys:
            tot[x]+=d[x]
            if d["region"]=="London": totlon[x]+=d[x]
      row=[daytodate(datetoday(date)-1)]# Go back a day because Zoe values are reported (and timestamped) the day after they occur
      for loc in locs:
        if loc=="London": src=totlon
        elif loc=="UK": src=tot
        else: src=dd[loc]
        row.append("%.4g"%(src["corrected_covid_positive"]/src["population"]*1e5))
      writer.writerow(row)
  selprev=pd.read_csv(selectedfn)
  graph=selprev.plot('Date',locs)
  graph.set_xticklabels(selprev.Date)
  #pyplot.rcParams["date.autoformatter.day"]="%Y-%m-%d"
  pyplot.show()

#if __name__=="__main__":
#  processdata("zoemapdata")
  
