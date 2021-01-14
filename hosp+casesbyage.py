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

url='https://api.coronavirus.data.gov.uk/v1/data?'
req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"cumAdmissionsByAge"}'; hospdata=get_data(url+req)['data']
req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"maleCases"}';          maledata=get_data(url+req)['data']
req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"femaleCases"}';        femaledata=get_data(url+req)['data']

def getdiff(data):
  newdata=[]
  day=datetoday(data[0]['date'])
  n=1
  while n<len(data) and datetoday(data[n]['date'])==day-n: n+=1# Find maximal contiguous date range
  for i in range(n-1,-1,-1):
    l={}
    for d in data[i]['blah']: l[d['age']]=d['value']
    if i<n-1:
      newage={age: l[age]-pl[age] for age in l}
      newage['date']=data[i]['date']
      newdata.append(newage)
    pl=l
  return newdata

newhosp=getdiff(hospdata)
newmale=getdiff(maledata)
newfemale=getdiff(femaledata)
assert(len(newmale)==len(newfemale) and newmale[0]['date']==newfemale[0]['date'])
newcases=[]
for (d,e) in zip(newmale,newfemale):
  f={}
  for age in d: f[age]=d[age]+e[age]
  f['date']=d['date']
  newcases.append(f)
newcases=newcases[:-1]# Last entry seems unreliable, I think because it using specimin date and there are biases with recent entries

def smooth(data):
  ages=[x for x in data[0].keys() if x!='date']
  n=len(data)
  smoothed=[]
  for i in range(n):
    d={'date': data[i]['date']}
    r=min(3,i,n-1-i)
    for age in ages:
      d[age]=sum(data[j][age] for j in range(i-r,i+r+1))/(2*r+1)
    smoothed.append(d)
  return smoothed

hosp=smooth(newhosp)
cases=smooth(newcases)
  
def makegraph(title='A graph', data=[], mindate='0000-00-00', ylabel='', outfn='temp.png'):
  po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
  
  # Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
  def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

  write('set terminal pngcairo font "sans,13" size 1920,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%outfn)
  write('set for [i=9:16] linetype i dashtype (20,7)')
  write('set key right')
  write('set title "%s"'%title)
  write('set ylabel "'+ylabel+'"')
  write('set xdata time')
  write('set format x "%Y-%m-%d"')
  write('set timefmt "%Y-%m-%d"')
  write('set tics scale 2,0.5')
  write('set xtics "2020-01-06", 604800')#%startdate)# Date labels on Mondays
  write('set xtics rotate by 45 right offset 0.5,0')
  write('set grid xtics ytics lc rgb "#dddddd" lt 1')
  write('set xtics nomirror')
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
  p.close();po.wait()
  print("Written graph to %s"%outfn)

def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  return tuple(int(y) for y in x.split("_to_"))
  
date=max(hosp[-1]['date'],cases[-1]['date'])
mindate=daytodate(datetoday(date)-90)
hospages=sorted((x for x in hosp[0] if x!='date'),key=lambda x:parseage(x)[0])
caseages=sorted((x for x in cases[0] if x!='date'),key=lambda x:parseage(x)[0])

data=[]
for age in ['18_to_64', '65_to_84', '85+']:
  data.append({
    'title': age.replace('_',' '),
    'values': [(d['date'],d[age]) for d in hosp]
  })
title='Hospital admissions for Covid-19 in England by age group. Source: https://coronavirus.data.gov.uk/ at '+date
makegraph(title=title, data=data, mindate=mindate, ylabel='Number of age group admitted', outfn='hospitaladmissionsbyage-abs.png')

title='Hospital admissions and confirmed cases ratios for Covid-19 in England. Source: https://coronavirus.data.gov.uk/ at '+date
cutoff0=65
cutoff1=80
lowages=[age for age in caseages if parseage(age)[0]>=cutoff0 and parseage(age)[1]<cutoff1]
highages=[age for age in caseages if parseage(age)[0]>=cutoff1]
data=[]
#data.append({
#  'title': 'Hospital admissions: #(aged 65+) / #(all ages)',
#  'values': [(d['date'],(d['65_to_84']+0*d['85+'])/sum(d[a] for a in hospages)*100) for d in hosp]
#})
data.append({
  'title': 'Hospital admissions: #(aged 85+) / #(aged 18-64 or 85+)',
  'values': [(d['date'],d['85+']/(d['18_to_64']+d['85+'])*100) for d in hosp]
})
data.append({
  'title': 'Confirmed cases: #(aged %d+) / #(aged %d+)'%(cutoff1,cutoff0),
  'values': [(d['date'],sum(d[a] for a in highages)/sum(d[a] for a in lowages+highages)*100) for d in cases if d['date']>=mindate]
})
if 0:
  num=[a for a in caseages if parseage(a)[0]>=85]
  denom=[a for a in caseages if parseage(a)[0]>=20 and parseage(a)[1]<=64 or parseage(a)[0]>=85]
  data.append({
    'title': 'Confirmed cases: #(aged 85+) / #(aged 20-64 or 85+)',
    'values': [(d['date'],sum(d[a] for a in num)/sum(d[a] for a in denom)*100) for d in cases if d['date']>=mindate]
  })
makegraph(title=title, data=data, mindate=mindate, ylabel='Percentage', outfn='admissionandcaseageratios.png')
