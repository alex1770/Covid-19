import time,calendar
from requests import get
from subprocess import Popen,PIPE

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def get_data(url):
  response = get(url, timeout=10)
  if response.status_code >= 400:
    raise RuntimeError(f'Request failed: { response.text }')
  return response.json()

ages=['0_to_5', '6_to_17', '18_to_64', '65_to_84', '85+']
url='https://api.coronavirus.data.gov.uk/v1/data?'
req='filters=areaType=nation;areaName=england&structure={"date":"date","daily":"cumAdmissionsByAge"}'
data=get_data(url+req)['data']
newdata=[]
n=len(data)
for i in range(n-1,-1,-1):
  l={}
  for d in data[i]['daily']:
    l[d['age']]=d['value']
  if i<n-1:
    newage={age: l[age]-pl[age] for age in l}
    newage['date']=data[i]['date']
    newdata.append(newage)
  pl=l

n=len(newdata)
smoothed=[]
for i in range(n):
  d={'date': newdata[i]['date']}
  r=min(3,i,n-1-i)
  for age in ages:
    d[age]=sum(newdata[j][age] for j in range(i-r,i+r+1))/(2*r+1)
  smoothed.append(d)
  
if 0:
  for d in smoothed:
    print(d['date'],end='')
    for age in ages: print(" %7.1f"%(d[age]),end='')
    print()

def makegraph(title='A graph', data=[], mindate='0000-00-00', ylabel='', outfn='temp.png'):
  p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
  
  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  write('set terminal pngcairo font "sans,13" size 1920,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%outfn)
  #write('set for [i=9:16] linetype i dashtype (20,7)')
  write('set key left')
  write('set title "%s"'%title)
  #write('set xlabel "Days since '+desc+perstring+' reached %g'%thr)
  write('set ylabel "'+ylabel+'"')
  write('set xdata time')
  write('set format x "%Y-%m-%d"')
  write('set timefmt "%Y-%m-%d"')
  write('set tics scale 2,0.5')
  write('set xtics "2020-01-06", 604800')#%startdate)# Date labels on Mondays
  write('set xtics rotate by 45 right offset 0.5,0')
  write('set grid xtics ytics lc rgb "#dddddd" lt 1')
  s='plot '
  first=True
  for dat in data:
    if not first: s+=', '
    first=False
    s+='"-" using 1:2 with lines lw 3 title "%s"'%(dat['title'])
  write(s)
  for dat in data:
    for (date,val) in dat['values']:
      if date>=mindate: write(date,val)
    write("e")
  p.close()
  print("Written graph to %s"%outfn)

date=smoothed[-1]['date']
mindate=daytodate(datetoday(date)-90)

data=[]
for age in ['18_to_64', '65_to_84', '85+']:
  data.append({
    'title': age.replace('_',' '),
    'values': [(d['date'],d[age]) for d in smoothed]
  })
title='Hospital admissions for Covid-19 in England by age group. Source: https://coronavirus.data.gov.uk/ at '+date
makegraph(title=title, data=data, mindate=mindate, ylabel='Number of age group admitted', outfn='hospitaladmissionsbyage-abs.png')

title='Hospital admissions for Covid-19 in England - proportions by age group. Source: https://coronavirus.data.gov.uk/ at '+date
data=[]
for age in ['18_to_64', '65_to_84', '85+']:
  data.append({
    'title': age.replace('_',' ')+' out of all ages',
    'values': [(d['date'],d[age]/sum(d[a] for a in ages)*100) for d in smoothed]
  })
data.append({
  'title': '85+ out of 18-64 and 85+',
  'values': [(d['date'],d['85+']/(d['18_to_64']+d['85+'])*100) for d in smoothed]
})
makegraph(title=title, data=data, mindate=mindate, ylabel='Percentage admitted', outfn='hospitaladmissionsbyage-rel.png')
