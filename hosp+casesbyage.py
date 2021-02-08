import time,calendar,os,json,sys,datetime
from requests import get
from subprocess import Popen,PIPE

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v1/data?'
  response = get(url+req, timeout=10)
  if not response.ok:
    raise RuntimeError(f'Request failed: { response.text }')
  date=time.strftime('%Y-%m-%d',time.strptime(response.headers['Last-Modified'],'%a, %d %b %Y %H:%M:%S %Z'))# Not currently used
  data=response.json()['data']
  
  # Convert from list form to dictionary keyed by age
  day=datetoday(data[0]['date'])
  n=1
  while n<len(data) and datetoday(data[n]['date'])==day-n: n+=1# Find maximal contiguous date range
  data1=[]
  for i in range(n-1,-1,-1):
    d=data[i]
    e={'date':d['date']}
    for x in d:
      if x!='date':
        for y in d[x]:
          e[y['age']]=e.get(y['age'],0)+y['value']
    data1.append(e)

  return data1

req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"cumAdmissionsByAge"}';               hospdata=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","male":"maleCases","female":"femaleCases"}'; casedata=get_data(req)
updatedate=casedata[-1]['date']
now=datetime.datetime.utcnow().strftime('%Y-%m-%d')

# Save case data because we might want to artificially implement cases-by-publication-date-and-age. (newCasesByPublishDateAgeDemographics not working)
fn=os.path.join('apidata',updatedate)
if len(sys.argv)==1 and os.path.isfile(fn): sys.exit(1)# Exit signalling no update needs to be done
os.makedirs('apidata', exist_ok=True)
with open(fn,'w') as fp:
  json.dump(casedata,fp,indent=2)

def getdiff(data):
  n=len(data)
  newdata=[]
  for i in range(1,n):
    l={'date':data[i]['date']}
    for age in data[i]:
      if age!='date': l[age]=data[i][age]-data[i-1][age]
    newdata.append(l)
  return newdata

newhosp=getdiff(hospdata)
newcases=getdiff(casedata)
newcases=newcases[:-1]# Last entry seems particularly unreliable, I think because it using specimen date and there are biases with recent entries

def smooth(data):
  ages=[x for x in data[0].keys() if x!='date']
  n=len(data)
  smoothed=[]
  for i in range(n):
    d={'date': data[i]['date']}
    j0=max(i-3,0)
    j1=min(i+4,n)
    for age in ages:
      d[age]=sum(data[j][age] for j in range(j0,j1))/(j1-j0)
    smoothed.append(d)
  return smoothed

hosp=smooth(newhosp)
cases=smooth(newcases)

def makegraph(title='A graph', data=[], mindate='0000-00-00', ylabel='', outfn='temp.png', extra=[]):
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
  for x in extra: write(x)
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
  
#date=max(hosp[-1]['date'],cases[-1]['date'])
mindate=daytodate(datetoday(updatedate)-90)
hospages=sorted((x for x in hosp[0] if x!='date'),key=lambda x:parseage(x)[0])
caseages=sorted((x for x in cases[0] if x!='date'),key=lambda x:parseage(x)[0])

data=[]
for age in ['18_to_64', '65_to_84', '85+']:
  data.append({
    'title': age.replace('_',' '),
    'values': [(d['date'],d[age]) for d in hosp]
  })
title='Hospital admissions for Covid-19 in England by age group. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
makegraph(title=title, data=data, mindate=mindate, ylabel='Number of age group admitted', outfn='hospitaladmissionsbyage-abs.png')

data=[]
for ageband in range(0,90,10):
  if ageband<80: lim=ageband+10;name="%d - %d"%(ageband,lim)
  else: lim=999;name="%d+"%ageband
  data.append({
    'title': name,
    'values': [(d['date'],sum(d[age] for age in caseages if parseage(age)[0]>=ageband and parseage(age)[1]<=lim)) for d in cases]
  })
title='Confirmed cases per day for Covid-19 in England by age group. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
makegraph(title=title, data=data, mindate=mindate, ylabel='Number of cases per day', outfn='confirmedcasesbyage-abs.png')#, extra=['set logscale y'])

title='Hospital admissions and confirmed cases ratios for Covid-19 in England. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
cutoff0=65
cutoff1=80
lowages=[age for age in caseages if parseage(age)[0]>=cutoff0 and parseage(age)[1]<cutoff1]
highages=[age for age in caseages if parseage(age)[0]>=cutoff1]
data=[]

if 0:
  data.append({
    'title': 'Hospital admissions: #(aged 65+) / #(all ages)',
    'values': [(d['date'],(d['65_to_84']+0*d['85+'])/sum(d[a] for a in hospages)*100) for d in hosp]
  })
  
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

