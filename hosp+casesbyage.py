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
          if 'value' in y: val=y['value']
          else: val=y['deaths']
          e[y['age']]=e.get(y['age'],0)+val
    data1.append(e)

  return data1

req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"newDeaths28DaysByDeathDateAgeDemographics"}'; mortdata=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","blah":"cumAdmissionsByAge"}';                        hospdata=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","male":"maleCases","female":"femaleCases"}';          casedata=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","male":"maleCases"}';              malecases=get_data(req)
req='filters=areaType=nation;areaName=england&structure={"date":"date","female":"femaleCases"}';          femalecases=get_data(req)
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
newmcases=getdiff(malecases)
newfcases=getdiff(femalecases)
newcases=newcases[:-1]# Last entry seems particularly unreliable, I think because it using specimen date and there are biases with recent entries
newmcases=newmcases[:-1]
newfcases=newfcases[:-1]
                  
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
deaths=smooth(mortdata)
mcases=smooth(newmcases)
fcases=smooth(newfcases)

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
    s+='"-" using 1:2 with lines '+dat.get('extra','')+' lw 3 title "%s"'%(dat['title'])
  write(s)
  for dat in data:
    for (date,val) in dat['values']:
      if date>=mindate: write(date,val)
    write("e")
  p.close();po.wait()
  print("Written graph to %s"%outfn)

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

# Convert (eg) (15,20) to "15 - 19"
def unparse(r):
  (a,b)=r
  if b==150: return "%d+"%a
  return "%d - %d"%(a,b)

# Convert dictionary from using '15_19' (etc) format to (15,20) format
# At the same time remove age brackets such as '60+' and '00_59' that strictly contain other age brackets, so avoiding overcounting
# Return list of ages
def convertages(dd):
  ages0=[(x,parseage(x)) for x in dd[0] if x!='date']
  ages1={}
  for (x,(a,b)) in ages0:
    for (y,(c,d)) in ages0:
      if c>=a and d<=b and (c>a or d<b): break
    else: ages1[x]=(a,b)
  ee=[]
  for d in dd:
    e={}
    e['date']=d['date']
    for x in ages1:
      e[ages1[x]]=d[x]
    ee.append(e)
  ages2=sorted(ages1.values())
  return (ee,ages2)

#date=max(hosp[-1]['date'],cases[-1]['date'])
mindate=daytodate(datetoday(updatedate)-90)
hosp,hospages=convertages(hosp)
cases,caseages=convertages(cases)
deaths,deathages=convertages(deaths)
fcases,_=convertages(fcases)
mcases,_=convertages(mcases)

data=[]
for age in [(18,65), (65,85), (85,150)]:
  data.append({
    'title': unparse(age),
    'values': [(d['date'],d[age]) for d in hosp]
  })
title='Hospital admissions for Covid-19 in England by age group. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
makegraph(title=title, data=data, mindate=mindate, ylabel='Number of age group admitted', outfn='hospitaladmissionsbyage-abs.png')

# Todo when can be bothered: normalise this by number in each age group
data=[]
for ageband in range(0,90,10):
  if ageband<80: lim=ageband+10
  else: lim=150
  data.append({
    'title': unparse((ageband,lim)),
    'values': [(d['date'],sum(d[age] for age in caseages if age[0]>=ageband and age[1]<=lim)) for d in cases]
  })
title='Confirmed cases per day for Covid-19 in England by age group. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
makegraph(title=title, data=data, mindate=mindate, ylabel='Number of cases per day', outfn='confirmedcasesbyage-abs.png')#, extra=['set logscale y'])

title='Hospital admissions and confirmed cases ratios for Covid-19 in England. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
cutoff0=65#25#65
cutoff1=80#75#80
lowages=[age for age in caseages if age[0]>=cutoff0 and age[1]<=cutoff1]
highages=[age for age in caseages if age[0]>=cutoff1]
data=[]

if 0:
  days=(range(330,340),[-1])
  l=[]
  #dd={}
  #for end in [0,1]:
  #  for cut in caseages+[150]:
  #    dd[(end,cut)]=sum(cases[day][age] for day in days[end] for age in caseages if parseage(age)[
  day0=336;day1=-1
  for cut0 in range(20,85,5):
    for cut1 in range(cut0+5,95,5):
      for cut2 in range(cut1,95,5):
        lowa=[age for age in caseages if age[0]>=cut0 and age[1]<=cut1]
        higha=[age for age in caseages if age[0]>=cut2]
        rr=[[sum(d[a] for a in aa) for aa in [lowa,higha]] for d in [cases[day0],cases[day1]]]
        if min(min(rr[0]),min(rr[1]))>=100:
          l.append((rr[0][1]/rr[0][0]/(rr[1][1]/rr[1][0]),cut0,cut1,cut2))
          #if cut0==65 and cut1==80: print(rr[0],rr[1])
  l.sort(reverse=True)
  for (r,cut0,cut1,cut2) in l[:20]:
    print(cut0,cut1,cut2,r)
  (r,cut0,cut1,cut2)=l[0]
  lowages=[age for age in caseages if age[0]>=cut0 and age[1]<cut1]
  highages=[age for age in caseages if age[0]>=cut2]

if 0:
  data.append({
    'title': 'Hospital admissions: #(aged 65+) / #(all ages)',
    'values': [(d['date'],(d[(65,85)]+0*d[(85,150)])/sum(d[a] for a in hospages)*100) for d in hosp]
  })
  
data.append({
  'title': 'Hospital admissions: #(aged 85+) / #(aged 18-64 or 85+)',
  'values': [(d['date'],(d[(85,150)])/(d[(18,65)]+d[(85,150)])*100) for d in hosp if d['date']>=mindate]
})

data.append({
  'title': 'Confirmed cases: #(aged %d+) / #(aged %d+)'%(cutoff1,cutoff0),
  'values': [(d['date'],sum(d[a] for a in highages)/sum(d[a] for a in lowages+highages)*100) for d in cases if d['date']>=mindate]
})

lowages=[age for age in deathages if age[0]>=cutoff0 and age[1]<=cutoff1]
highages=[age for age in deathages if age[0]>=cutoff1]
data.append({
  'title': 'Deaths: #(aged %d+) / #(aged %d+) - 25%%'%(cutoff1,cutoff0),
  'values': [(d['date'],sum(d[a] for a in highages)/sum(d[a] for a in lowages+highages)*100-25) for d in deaths if d['date']>=mindate],
  #'extra': 'axis x1y2'
})

if 0:
  num=[a for a in caseages if a[0]>=85]
  denom=[a for a in caseages if a[0]>=20 and a[1]<=65 or a[0]>=85]
  data.append({
    'title': 'Confirmed cases: #(aged 85+) / #(aged 20-64 or 85+)',
    'values': [(d['date'],sum(d[a] for a in num)/sum(d[a] for a in denom)*100) for d in cases if d['date']>=mindate]
  })

makegraph(title=title, data=data, mindate=mindate, ylabel='Percentage', outfn='admissionandcaseageratios.png')

data=[]
lowages=[age for age in caseages if age[0]>=16 and age[1]<=65]
data.append({
  'title': 'Confirmed cases: #(female aged 16-65) / #(male aged 16-65)',
  'values': [(f['date'],sum(f[a] for a in lowages)/sum(m[a] for a in lowages)) for (f,m) in zip(fcases,mcases) if f['date']>=mindate]
})
makegraph(title=title, data=data, mindate=mindate, ylabel='Ratio', outfn='femalemalecaseratio.png')
