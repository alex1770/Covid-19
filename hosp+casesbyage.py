import time,calendar,os,json,sys,datetime
from requests import get
from subprocess import Popen,PIPE
from math import sqrt,log
from scipy.optimize import minimize
import numpy as np

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
                  
def simplesmooth(data):
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
def smooth(data):
  return simplesmooth(data)

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
hosp,hospages=convertages(newhosp)
cases,caseages=convertages(newcases)
deaths,deathages=convertages(mortdata)
fcases,_=convertages(newfcases)
mcases,_=convertages(newmcases)

def LL(rr,xx):
  L=0
  n=len(rr)
  for i in range(7):
    x=xx[i::7].sum()
    w=x/(rr[i::7].sum())
    L+=x*log(w)
  L+=(xx*np.log(rr)).sum()
  dd=-rr[:-2]+2*rr[1:-1]-rr[2:]
  t=(dd*dd).sum()
  s=(rr*rr).sum()
  #print(s,t)
  #L-=n*log(t/s)
  lam=100000000;L-=lam*t
  return -L/100000

def val2(rr,xx,vv,lam):
  dd=-rr[:-2]+2*rr[1:-1]-rr[2:]
  return -(vv*rr).sum()+(xx*np.log(vv*rr)).sum()-lam/2*(dd*dd).sum()

def magicmethod(xx):
  n=len(xx)
  ww=[xx[i::7].sum()/len(xx[i::7]) for i in range(7)]
  rr=np.array([xx[d]/ww[d%7] for d in range(n)])
  rr=rr/rr.sum()

  lam=1
  
  while 1:
    ww=[xx[i::7].sum()/rr[i::7].sum() for i in range(7)]
    vv=[ww[d%7] for d in range(n)]
    print(val2(rr,xx,vv,lam))

    A=np.zeros((n,n))
    for d in range(1,n-1):
      A[d-1,d-1]+=1
      A[d-1,d  ]-=2
      A[d-1,d+1]+=1
      A[d  ,d-1]-=2
      A[d  ,d  ]+=4
      A[d  ,d+1]-=2
      A[d+1,d-1]+=1
      A[d+1,d  ]-=2
      A[d+1,d+1]+=1
    A=lam/2*A
    B=vv-xx/rr
    res=np.linalg.lstsq(A,B)
    rr=res[0]
    a=sqrt((rr*rr).sum()/n)
    rr=np.maximum(rr,a/100)
    
    
if 1:
  data=cases
  ages=[x for x in data[0].keys() if x!='date']
  xx=np.array([sum(d[age] for age in ages) for d in data])
  n=len(data)
  ww=[xx[i::7].sum()/len(xx[i::7]) for i in range(7)]
  rr=np.array([xx[d]/ww[d%7] for d in range(n)])
  rr=rr/sqrt((rr*rr).sum())
  constr={'type':'eq', 'fun':lambda rr: (rr*rr).sum()-1}
  res=minimize(LL,rr,args=(xx),method="SLSQP",bounds=[(1e-9,None) for i in range(n)],constraints=[constr],options={"maxiter":1000})
  print(res.success)
  print(res.message)
  print(res.nit,"iterations")
  rr=res.x
  print(LL(rr,xx))
  print()

  n=len(xx)
  ww=[xx[i::7].sum()/rr[i::7].sum() for i in range(7)]
  vv=[ww[d%7] for d in range(n)]
  print((-(vv*rr).sum()))
  print((xx*np.log(vv*rr)).sum())
  dd=-rr[:-2]+2*rr[1:-1]-rr[2:]
  t=(dd*dd).sum()
  s=(rr*rr).sum()
  print(t,s,n*log(t/s))
  with open('temp','w') as fp:
    for i in range(n):
      print("%12g %12g %12g"%(xx[i],rr[i],vv[i]),file=fp)
  
  #rr=magicmethod(xx)
  
  smoothed=[]
  for i in range(n):
    d={'date': data[i]['date']}
    j0=max(i-3,0)
    j1=min(i+4,n)
    for age in ages:
      d[age]=sum(data[j][age] for j in range(j0,j1))/(j1-j0)
    smoothed.append(d)
  #return smoothed
  poiopi


hosp=smooth(hosp)
cases=smooth(cases)
deaths=smooth(deaths)
mcases=smooth(mcases)
fcases=smooth(fcases)

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

if 0:
  days=(range(330,340),[-1])
  ll=[]
  for (ages,numthings,desc) in [(caseages,cases,"cases"), (deathages,deaths,"deaths")]:
    aa={}
    dd={}
    for end in [0,1]:
      for cut in [x[0] for x in ages]+[150]:
        dd[(end,cut)]=sum(numthings[day][age] for day in days[end] for age in ages if age[0]<cut)/len(days[end])
    n=len(ages)
    for c0 in range(n-2):
      cut0=ages[c0][0]
      for c1 in range(c0+1,n-1):
        cut1=ages[c1][0]
        for c2 in range(c1,n):
          cut2=ages[c2][0]
          rr=[]
          for end in [0,1]:
            rr.append(dd[(end,cut1)]-dd[(end,cut0)])
            rr.append(dd[(end,150)] -dd[(end,cut2)])
          if min(rr)>=10:
            aa[cut0,cut1,cut2]=rr[1]/rr[0]/(rr[3]/rr[2])
    ll.append(aa)
  l=[]
  for x in ll[0]:
    if x in ll[1]:
      l.append((sqrt(ll[0][x]*ll[1][x]),*x))
  l.sort(reverse=True)
  for (r,cut0,cut1,cut2) in l:
    if cut2<=70: print("%2d %2d %2d %7.3f"%(cut0,cut1,cut2,r))
    if r<0.9*l[0][0]: break


title='Hospital admissions and confirmed cases/deaths ratios for Covid-19 in England. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
data=[]

cutoff0=0;cutoff1=65;cutoff2=65
lowages=[age for age in hospages if age[0]>=cutoff0 and age[1]<=cutoff1]
highages=[age for age in hospages if age[0]>=cutoff2]
data.append({
  'title': 'Hospital admissions: 0.75 * (aged %s) / (aged %s)'%(unparse((highages[0][0],highages[-1][1])),unparse((lowages[0][0],lowages[-1][1]))),
  'values': [(d['date'],0.75*sum(d[a] for a in highages)/sum(d[a] for a in lowages)) for d in hosp if d['date']>=mindate]
})

cutoff0=0;cutoff1=60;cutoff2=70
lowages=[age for age in caseages if age[0]>=cutoff0 and age[1]<=cutoff1]
highages=[age for age in caseages if age[0]>=cutoff2]
data.append({
  'title': 'Confirmed cases: 10 * (aged %s) / (aged %s)'%(unparse((highages[0][0],highages[-1][1])),unparse((lowages[0][0],lowages[-1][1]))),
  'values': [(d['date'],10*sum(d[a] for a in highages)/sum(d[a] for a in lowages)) for d in cases if d['date']>=mindate]
})

lowages=[age for age in deathages if age[0]>=cutoff0 and age[1]<=cutoff1]
highages=[age for age in deathages if age[0]>=cutoff2]
data.append({
  'title': 'Deaths: 0.1 * (aged %s) / (aged %s)'%(unparse((highages[0][0],highages[-1][1])),unparse((lowages[0][0],lowages[-1][1]))),
  'values': [(d['date'],0.1*sum(d[a] for a in highages)/sum(d[a] for a in lowages)) for d in deaths if d['date']>=mindate],
  #'extra': 'axis x1y2'
})

makegraph(title=title, data=data, mindate=mindate, ylabel='Adjusted Ratio', outfn='admissionandcaseageratios2.png')


#################################
# Old graphs (14 Jan - 5 March) #
#################################
    
title='Hospital admissions and confirmed cases/deaths ratios for Covid-19 in England. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
cutoff0=65;cutoff1=150;cutoff2=80
data=[]

data.append({
  'title': 'Hospital admissions: (aged 85+) / (aged 18-64 or 85+)',
  'values': [(d['date'],(d[(85,150)])/(d[(18,65)]+d[(85,150)])*100) for d in hosp if d['date']>=mindate]
})

lowages=[age for age in caseages if age[0]>=cutoff0 and age[1]<=cutoff1]
highages=[age for age in caseages if age[0]>=cutoff2]
data.append({
  'title': 'Confirmed cases: (aged %s) / (aged %s)'%(unparse((cutoff2,150)),unparse((cutoff0,cutoff1))),
  'values': [(d['date'],sum(d[a] for a in highages)/sum(d[a] for a in lowages)*100) for d in cases if d['date']>=mindate]
})

lowages=[age for age in deathages if age[0]>=cutoff0 and age[1]<=cutoff1]
highages=[age for age in deathages if age[0]>=cutoff2]
data.append({
  'title': 'Deaths: (aged %s) / (aged %s) - 25%%'%(unparse((cutoff2,150)),unparse((cutoff0,cutoff1))),
  'values': [(d['date'],sum(d[a] for a in highages)/sum(d[a] for a in lowages)*100-25) for d in deaths if d['date']>=mindate],
  #'extra': 'axis x1y2'
})

makegraph(title=title, data=data, mindate=mindate, ylabel='Percentage', outfn='admissionandcaseageratios.png')

########################

data=[]
lowages=[age for age in caseages if age[0]>=16 and age[1]<=65]
data.append({
  'title': 'Confirmed cases: #(female aged 16-65) / #(male aged 16-65)',
  'values': [(f['date'],sum(f[a] for a in lowages)/sum(m[a] for a in lowages)) for (f,m) in zip(fcases,mcases) if f['date']>=mindate]
})
makegraph(title=title, data=data, mindate=mindate, ylabel='Ratio', outfn='femalemalecaseratio.png')

########################

data=[]
for age in [(18,65), (65,85), (85,150)]:
  data.append({
    'title': unparse(age),
    'values': [(d['date'],d[age]) for d in hosp]
  })
title='Hospital admissions for Covid-19 in England by age group. Last few values subject to change.\\nSource: https://coronavirus.data.gov.uk/ at '+now
makegraph(title=title, data=data, mindate=mindate, ylabel='Number of age group admitted', outfn='hospitaladmissionsbyage-abs.png')

########################

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

