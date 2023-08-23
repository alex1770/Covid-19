import csv,time,calendar,os,json,sys,datetime,subprocess

# Input: 'it' is an iterable (such as a file pointer)
# Returns map: {headings} -> [list of entries]
def loadcsv_it(it,sep=","):
  reader=csv.reader(it,delimiter=sep)
  headings=next(reader)
  out={}
  for (i,h) in enumerate(headings):
    if h in out: headings[i]=None
    else: out[h]=[]
  for row in reader:
    #assert len(row)==len(headings)
    for (h,x) in zip(headings,row):
      if h!=None: out[h].append(x)
  for h in out:
    try:
      out[h]=[int(x.replace(',','')) for x in out[h]]
    except ValueError:
      try:
        out[h]=[float(x.replace(',','')) for x in out[h]]
      except ValueError:
        pass
  return out

def loadcsv(fn,sep=","):
  with open(fn,'r') as fp:
    return loadcsv_it(fp,sep)

# Generator for rows of csv with specified headings
# (No type conversion - all items will be returned as strings)
# If the heading is None then the entire line will be returned
def csvrows_it(it,reqheadings,sep=","):
  reader=csv.reader(it,delimiter=sep)
  headings=next(reader)
  dec={None:None}
  for (i,h) in enumerate(headings): dec[h]=i
  cols=[]
  for h in reqheadings:
    if h in dec:
      cols.append(dec[h])
    else:
      raise RuntimeError('Heading '+h+' not found in csv file')
  while 1:
    row=next(reader,None)
    if row==None: return
    yield [(row if i==None else row[i]) for i in cols]
  
def csvrows(fn,reqheadings,sep=','):
  with open(fn,'r') as fp:
    for x in csvrows_it(fp,reqheadings,sep):
      yield x

def savecsv(incsv,fn,reqheadings=None,sep=','):
  if reqheadings==None: reqheadings=list(incsv)
  with open(fn,'w') as fp:
    writer=csv.writer(fp)
    writer.writerow(list(incsv))
    for row in zip(*[incsv[h] for h in reqheadings]):
      writer.writerow(row)

def datetoday(x):
  format=None
  if x[-2:]=='US':
    format='%m/%d/%YUS';# E.g., 06/25/2021US
  elif '/' in x:
    format='%d/%m/%Y';# E.g., 25/06/2021
  elif ' ' in x:
    if len(x)==11:
      format='%d %b %Y'# E.g., 05 Jun 2021
    else:
      format='%d %B %Y'# E.g., 05 June 2021
  elif '-' in x:
    y=x.split('-')
    if len(y[0])==4 and y[1].isdigit():
      format='%Y-%m-%d';# E.g., 2021-06-05
    elif y[1].isalpha() and len(y[2])==2:
      format='%d-%b-%y';# E.g., 29-Oct-20
  elif len(x)==8 and x.isdigit():
    format="%Y%m%d"
  if format==None: raise RuntimeError('Unrecognised date format: '+x)
  t=time.strptime(x+'UTC',format+'%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

def api_v2(req):
  import requests
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  response=requests.get(url+req,timeout=10)
  if not response.ok:
    raise RuntimeError('Request failed: '+response.text)
  return response

def UKdatetime():
  import pytz
  now=datetime.datetime.now(tz=pytz.timezone('Europe/London'))
  return Date(now.strftime('%Y-%m-%d')),now

def apiday():
  import pytz
  nowdate,nowtime=UKdatetime()
  # Dashboard is updated at 4pm UK time, so go back a day if time in UK < 4pm
  date=nowdate-int(nowtime.hour<16)
  
  # Intermediate case: return the final Friday update of 2022-07-01
  if date>="2022-07-01" and date<"2022-07-06": return Date("2022-07-01")
  
  # From 2022-07-04, dashboard is now only updated on Wednesdays, so after 2022-07-06, go back to previous Wednesday:
  w=(nowtime.weekday()-int(nowtime.hour<16)-2)%7# Weekday relative to Wednesday=0
  return date-w

def makegraph(title='A graph', data=[], mindate='0000-00-00', ylabel='', outfn='temp.png', extra=[], interval=604800, ranges=''):
  po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE);p=po.stdin
  
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
  write('set xtics "2020-01-06", %d'%interval)# Date labels on Mondays
  write('set xtics rotate by 45 right offset 0.5,0')
  write('set grid xtics ytics lc rgb "#dddddd" lt 1')
  write('set xtics nomirror')
  for x in extra: write(x)
  s='plot '+ranges+' '
  first=True
  for dat in data:
    if not first: s+=', '
    first=False
    (plotwith,n)=dat.get('with',('lines',1))
    s+='"-" using 1'
    for i in range(2,n+2): s+=':%d'%i
    s+=' with '+plotwith
    s+=' '+dat.get('extra','')
    s+=' lw 3 title "%s"'%(dat['title'])
  write(s)
  for dat in data:
    for dv in dat['values']:
      if dv[0]>=mindate and dv[1]!=None: write(*dv)
    write("e")
  p.close();po.wait()
  print("Written graph to %s"%outfn)

# Simple weekday adjustment
def weekdayadj(nn,eps=0.5):
  # y_{day} = log(nn(day)+eps) + w_{day%7}
  # error = sum_{day} (y_{day+1}-y_{day})^2
  # Choose w0, ..., w5 (w6 = 0, gauge fix) to minimise error
  from math import log,exp
  import numpy as np

  n=len(nn)
  targ=[log(x+eps) for x in nn]
  A=np.zeros([6,6])
  b=np.zeros(6)
  c=0
  for i in range(n-1):
    d=targ[i+1]-targ[i]
    i0=i%7
    i1=(i+1)%7
    if i1<6: A[i1,i1]+=1
    if i0<6: A[i0,i0]+=1
    if i0<6 and i1<6: A[i0,i1]-=1;A[i1,i0]-=1
    if i1<6: b[i1]-=d
    if i0<6: b[i0]+=d
    c+=d*d
  if np.__version__<'1.14':
    ww=np.linalg.lstsq(A,b)[0]
  else:
    ww=np.linalg.lstsq(A,b,rcond=None)[0]
  ww7=list(ww)+[0]
  weekadj=[exp(x) for x in ww7]
  adjusted=[nn[d]*weekadj[d%7] for d in range(n)]
  s=sum(nn)/sum(adjusted)
  return np.array([s*x for x in adjusted])

# Weekday adjustment incorportaing holidays
# Not ready - may be deleted
def weekdayadj2(nn,date0=None,eps=0.5):
  # y_{day} = log(nn(day)+eps) + w_{day%7}
  # error = sum_{day} (y_{day+1}-y_{day})^2
  # Choose w0, ..., w5 to minimise this error (w6 = 0, gauge fix)
  monday=Date("2022-01-03")
  holidays={
    "2020-01-01","2020-04-10","2020-04-13","2020-05-08","2020-05-25",             "2020-08-31","2020-12-25","2020-12-28",
    "2021-01-01","2021-04-02","2021-04-05","2021-05-03","2021-05-31",             "2021-08-30","2021-12-27","2021-12-28",
    "2022-01-03","2022-04-15","2022-04-18","2022-05-02","2022-06-02","2022-06-03","2022-08-29","2022-12-26","2022-12-27",
    "2023-01-02","2023-04-07","2023-04-10","2023-05-01","2023-05-29",             "2023-08-28","2023-12-25","2023-12-26"
  }
  from math import log,exp
  import numpy as np

  for holiday_day_of_week in range(8):
    if holiday_day_of_week==7: date0=None
    n=len(nn)
    targ=[log(x+eps) for x in nn]
    if date0 is None:
      day=[i%7 for i in range(n)]
    else:
      day=[]
      date0=Date(date0)
      for i in range(n):
        date=date0+i
        if str(date) in holidays: day.append(holiday_day_of_week)
        else: day.append((date-monday)%7)
    A=np.zeros([6,6])
    b=np.zeros(6)
    for i in range(n-1):
      diff=targ[i+1]-targ[i]
      i0=day[i]
      i1=day[i+1]
      if i1<6: A[i1,i1]+=1
      if i0<6: A[i0,i0]+=1
      if i0<6 and i1<6: A[i0,i1]-=1;A[i1,i0]-=1
      if i1<6: b[i1]-=diff
      if i0<6: b[i0]+=diff
    if np.__version__<'1.14':
      ww=np.linalg.lstsq(A,b)[0]
    else:
      ww=np.linalg.lstsq(A,b,rcond=None)[0]
    ww7=list(ww)+[0]
    err0=err1=0
    for i in range(n-1):
      err0+=(targ[i+1]-targ[i])**2
      err1+=(targ[i+1]+ww7[day[i+1]]-(targ[i]+ww7[day[i]]))**2
    print("Errors",holiday_day_of_week,err0,err1)
  weekadj=[exp(x) for x in ww7]
  adjusted=[nn[i]*weekadj[day[i]] for i in range(n)]
  s=sum(nn)/sum(adjusted)
  return np.array([s*x for x in adjusted])

# Slower weekday adjustment (n parameters)
def weekdayadj_slow(nn,alpha=0.1):
  # Find w0,...,w6 and Poisson parameters lambda_i to maximise
  # Sum_i log(P(Poisson(w_{i%7}*lambda_i)=nn_i)) - alpha*sum_i (log(lambda_{i+1})-log(lambda_i))^2
  from scipy.optimize import minimize
  from math import log
  import numpy as np
  n=len(nn)
  nn=np.array(nn)
  def LL(xx):
    ll=(nn*xx).sum()
    ee=np.exp(xx)
    for d in range(7):
      ll-=nn[d::7].sum()*log(ee[d::7].sum())
    dd=xx[1:]-xx[:-1]
    ll-=alpha*(dd*dd).sum()
    return ll
  xx=np.log(np.array(nn)+.5)
  LL0=LL(xx)# Baseline so that NLL returns small values
  def NLL(xx): return LL0-LL(xx)
  bounds=[(x-2,x+2) for x in xx]
  bounds[-1]=(xx[-1],xx[-1])# Gauge fix because NLL is invariant under xx -> xx+const
  res=minimize(NLL,xx,bounds=bounds,method="SLSQP",options={'maxiter':10000})#, 'eps':1e-4, 'ftol':1e-12})
  if not res.success: raise RuntimeError(res.message)
  adjusted=np.exp(res.x)
  return adjusted*nn.sum()/adjusted.sum()

# Not doing strict typechecking for all operations: e.g., -Date('2021-01-01') = -18628, Date('2021-01-01')+Date('2021-01-01') = Date('2072-01-02')
# Can't remember all the reasons, but I think the advantage of making it a subclass of int is that you don't have to define arithmetic operations like *, /, %, etc: it will just decay to an int.
class Date(int):
  def __new__(cls, daydate):
    if type(daydate)==str: return super(cls,cls).__new__(cls,datetoday(daydate))
    elif type(daydate)==datetime.datetime: return super(cls,cls).__new__(cls,datetoday(daydate.strftime("%Y-%m-%d")))
    else: return super(cls,cls).__new__(cls,daydate)
  def __eq__(self,other): return int(self)==int(Date(other))
  def __hash__(self): return hash(int(Date(self)))
  def __ge__(self,other): return int(self)>=int(Date(other))
  def __le__(self,other): return int(self)<=int(Date(other))
  def __gt__(self,other): return int(self)>int(Date(other))
  def __lt__(self,other): return int(self)<int(Date(other))
  # Date + int = Date (arithmetic on days)
  # Date + str = str (string concatenation)
  def __add__(self, other):
    if type(other)==str: return daytodate(int(self))+other
    # Can't remember why we don't just return Date(int(self)+int(other)) here
    res = super(Date, self).__add__(int(other))# int(other) rather than other makes it work when other is a numpy.int64
    return self.__class__(res)
  def __radd__(self, other):
    if type(other)==str: return other+daytodate(int(self))
    return self+other
  # Date - int = Date
  # Date - Date = int
  def __sub__(self, other):
    if type(other)==int: return Date(int(self)-other)
    return int(self)-int(Date(other))
  def __rsub__(self, other): return -(self-other)
  def __str__(self): return daytodate(int(self))
  def __repr__(self): return 'Date('+daytodate(int(self))+')'

def Daterange(a,b,step=1):
  a1=Date(a)
  b1=Date(b)
  for x in range(a1,b1,step):
    yield Date(x)

def gettopdir():
  f=__file__
  while 1:
    try:
      f=os.path.realpath(os.readlink(f))
    except OSError:
      break
  return os.path.dirname(f)

def getapidata(req,failsoft=False):
  import requests
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(5):
    try:
      response = requests.get(url+req, timeout=5)
      if response.ok: break
      error=response.text
    except BaseException as err:
      error=str(err)
  else:
    if failsoft: return None
    raise RuntimeError('Request failed: '+error)
  return response.json()['body'][::-1]

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('-','_').replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

def getpublishdate():
  t=getapidata('areaType=nation&areaName=England&metric=newCasesByPublishDate')
  return Date(t[-1]['date'])

def getcasesbyage_raw(pubday,location):
  # Target save format is
  # filename=publishdate, td[sex][specimendate][agerange] = cumulative cases,
  # (td:sex -> specdate -> agestring -> number_of_cases)
  # having converted agerange to open-closed format and eliminated superfluous ranges, but kept as a string because json can't handle tuples
  # Note that specimendate goes back to the dawn of time
  origages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
  cachedir=os.path.join(gettopdir(),'apidata_allcaseages')
  if location=='England':
    areatype='nation'
  else:
    areatype='region'
    cachedir+='_'+location
  date=str(Date(pubday))
  fn=os.path.join(cachedir,date)
  os.makedirs(cachedir,exist_ok=True)
  if os.path.isfile(fn):
    with open(fn,'r') as fp: td=json.load(fp)
  else:
    male=getapidata('areaType='+areatype+'&areaName='+location+'&metric=maleCases&release='+date)
    female=getapidata('areaType='+areatype+'&areaName='+location+'&metric=femaleCases&release='+date)
    td={}
    for sex in [male,female]:
      sexname=sex[0]['metric'][:-5]
      td[sexname]={}
      for d in sex:
        specdate=d['date']
        td[sexname][specdate]={}
        x=d[d['metric']]
        for y in x:
          a=parseage(y['age'])
          if a in origages:
            td[sexname][specdate]["%d_%d"%a]=y['value']
    with open(fn,'w') as fp: json.dump(td,fp,indent=2)
    print("Retrieved (fe)maleCases api data at",date,"for",location)
  return td

def getcases_raw(pubday,location="England"):
  # Target save format is
  # filename=publishdate, td[specimendate] = new cases,
  cachedir=os.path.join(gettopdir(),'apidata_allcases')
  if location=='England':
    areatype='nation'
  else:
    areatype='region'
    cachedir+='_'+location
  date=str(Date(pubday))
  fn=os.path.join(cachedir,date)
  os.makedirs(cachedir,exist_ok=True)
  if os.path.isfile(fn):
    with open(fn,'r') as fp: td=json.load(fp)
    if "Bad" in td: return td
  else:
    data=getapidata('areaType='+areatype+'&areaName='+location+'&metric=newCasesBySpecimenDate&release='+date,failsoft=True)
    if data==None:
      td={"Bad":True}
      print("Data not available from api at",date,": newCasesBySpecimenDate for",location)
      if pubday>=apiday(): return td# Don't permanently save this unavailability if it might become available later
    else:
      td={}
      for item in data:
        td[item["date"]]=int(item["newCasesBySpecimenDate"])
        td["Comment"]="newCasesBySpecimenDate"
      print("Retrieved newCasesBySpecimenDate api data published on",date,"for",location)
    with open(fn,'w') as fp: json.dump(td,fp,indent=2)
  if "Comment" in td: del td["Comment"]
  if "Bad" in td: return td
  return {Date(d):c for (d,c) in td.items()}

def getvirustests_raw(pubday,location="England"):
  def myint(x):
    try:
      return int(x)
    except:
      return None
  # Target save format is
  # filename=publishdate, td[specimendate] = tests
  cachedir=os.path.join(gettopdir(),'apidata_tests')
  if location=='England':
    areatype='nation'
  else:
    areatype='region'
    cachedir+='_'+location
  date=str(Date(pubday))
  fn=os.path.join(cachedir,date)
  os.makedirs(cachedir,exist_ok=True)
  if os.path.isfile(fn):
    with open(fn,'r') as fp: td=json.load(fp)
  else:
    data=getapidata('areaType='+areatype+'&areaName='+location+'&metric=newVirusTestsByPublishDate&metric=newVirusTestsBySpecimenDate&release='+date,failsoft=True)
    if data==None:
      td={"Bad":True}
      print("Data not available from api at",date,": newVirusTestsBy<Publish/Specimen>Date for",location)
      if pubday>=apiday(): return td# Don't permanently save this unavailability if it might become available later
    else:
      td={}
      for item in data:
        td[item["date"]]=[myint(item["newVirusTestsByPublishDate"]),myint(item["newVirusTestsBySpecimenDate"])]
        td["Comment"]=["newVirusTestsByPublishDate","newVirusTestsBySpecimenDate"]
      print("Retrieved newVirusTestsBy<Publish/Specimen>Date api data at",date,"for",location)
    with open(fn,'w') as fp: json.dump(td,fp,indent=2)
  if "Comment" in td: del td["Comment"]
  return td


# Returns numpy arrays:
#
# cc[pastpublishday - (minday-1)][specimenday - (minday-1)][sex 0=m, 1=f][index into ages] = pastpublishday's version of cumulative cases up to specimen day
# for minday-1 <= pastpublishday <= maxday, minday-1 <= specimenday <= maxday-1
#
# cn[pastpublishday - (minday-1)][specimenday - minday][sex 0=m, 1=f][index into ages] = pastpublishday's version of new cases on specimen day
# for minday-1 <= pastpublishday <= maxday, minday <= specimenday <= maxday-1
# 
# nn[pastpublishday - minday][specimenday - minday][sex 0=m, 1=f][index into ages] = new cases on specimen day that were first reported on pastpublish day
# for minday <= pastpublishday <= maxday, minday <= specimenday <= maxday-1
# (Some nn[] values can be <0 due to anomalous cumulative counts)
#
def convcasesbyagetonumpy(dd,minday,maxday,ages=[(0,150)]):
  import numpy as np
  origages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
  minday=Date(minday)
  maxday=Date(maxday)
  astrings=["%d_%d"%a for a in origages]
  reduceages={}
  for (a,astr) in enumerate(astrings):
    for (da,dat) in enumerate(ages):
      if origages[a][0]>=dat[0] and origages[a][1]<=dat[1]: reduceages[astr]=da;break
  nages=len(ages)
  
  # Convert to numpy array
  # cc[pastpublishday - (minday-1)][specimenday - (minday-1)][sex 0=m, 1=f][index into ages] = pastpublishday's version of cumulative cases up to specimen day
  # for minday-1 <= pastpublishday <= maxday, minday-1 <= specimenday <= maxday-1
  npub=maxday-minday+1
  nspec=maxday-minday
  cc=np.zeros([npub+1,nspec+1,2,nages],dtype=int)
  smindate=str(minday-1)# Prepare this to compare strings because datetoday is slow
  for pubdate in dd:
    pday=int(pubdate)-(minday-1)
    assert pday>=0
    for sex in dd[pubdate]:
      s=['male','female'].index(sex)
      for specdate in dd[pubdate][sex]:
        if specdate>=smindate:
          sday=specdate-(minday-1)
          assert sday>=0
          if sday<nspec+1:
            for astring in dd[pubdate][sex][specdate]:
              if astring in reduceages:
                cc[pday][sday][s][reduceages[astring]]+=dd[pubdate][sex][specdate][astring]
  
  # cn[pastpublishday - (minday-1)][specimenday - minday][sex 0=m, 1=f][index into ages] = pastpublishday's version of new cases on specimen day
  # for minday-1 <= pastpublishday <= maxday, minday <= specimenday <= maxday-1
  cn=cc[:,1:,:,:]-cc[:,:-1,:,:]
  for i in range(nspec): cn[i+1,i,:,:]=0

  # nn[pastpublishday - minday][specimenday - minday][sex 0=m, 1=f][index into ages] = new cases on specimen day that were first reported on pastpublish day
  # for minday <= pastpublishday <= maxday, minday <= specimenday <= maxday-1
  nn=cn[1:,:,:,:]-cn[:-1,:,:,:]
  
  return npub,nspec,cc,cn,nn

def getcasesbyagepubspec(minday,maxday,ages=[(0,150)],location='England'):
  # Get cumulative age-pubdate-specdate values from api or cache files
  dd={}
  for day in Daterange(minday-1,maxday+1):
    dd[day]=getcasesbyage_raw(day,location)

  # Convert to numpy array, taking difference of cumulative values to get incremental values
  return convcasesbyagetonumpy(dd,minday,maxday,ages=ages)

# Return case counts by age and specimen day, and estimated "completed" counts.
# Return values are
#   sp0[:,:], sp[:,:]
# where
#   sp0[day-minday,ageind] = number of cases for age index ageind on day 'day'
#    sp[day-minday,ageind] = estimated final number of cases for age index ageind on day 'day'
#   minday <= day < maxday
def getcasesbyagespeccomplete(maxday,ages=[(a,a+10) for a in range(0,70,10)]+[(70,150)],minday=Date('2021-08-20'),location='England'):
  infinity=7
  import numpy as np

  nages=len(ages)
  maxday=Date(maxday)
  
  dd={}
  for day in Daterange(maxday-max(7,infinity-2),maxday+1):
    dd[day]=getcasesbyage_raw(day,location)

  npub,nspec,cc,cn,nn=convcasesbyagetonumpy(dd,minday,maxday,ages=ages)
  gg=cn.sum(axis=2)# Sum over sexes
  
  # Try to undo the effect of delay from specimen to published test result by assuming the pattern is the same as last week's
  # and also adjusting for more recent trends in delayed reporting
  # sp[specimenday-minday][age index] = Est no. of samples
  sp=np.zeros([nspec,nages],dtype=float)
  for i in range(nspec):
    if npub-(i+1)>=infinity:
      sp[i]=gg[npub,i,:]
    else:
      base=np.array([gg[npub-r,i-r] for r in range(max(8,infinity-(npub-(i+1))))])
      targ7=gg[i+infinity-7,i-7]
      f0=1+sum((gg[npub,i-r]-gg[npub-1,i-r])/base[r] for r in range(1,infinity-(npub-(i+1))))
      f1=targ7/base[7]
      sp[i]=base[0]*(0.55*f0+0.45*f1)

  return gg[npub,:,:],sp

# Return (age-averaged) incomplete sample-day correction factors (between 0 and 1) in an array whose
# sample days correspond to (maxday-n, maxday-(n-1), ..., maxday-1).
def getextrap(maxday,ages=[(a,a+10) for a in range(0,70,10)]+[(70,150)],minday=Date('2021-08-20'),location='England'):
  sp0,sp=getcasesbyagespeccomplete(maxday,ages=ages,minday=minday,location=location)
  return sp0.sum(axis=1)/sp.sum(axis=1)

# https://gensplore.genomium.org/?gb=%2Fsequence.gb
# Counting from 1, inclusive-exclusive notation
genes={
  'ORF1a': [266, 13469],
  'ORF1b': [13468, 21556],
  'S': [21563, 25385],
  'ORF3a': [25393, 26221],
  'E': [26245, 26473],
  'M': [26523, 27192],
  'ORF6': [27202, 27388],
  'ORF7a': [27394, 27760],
  'ORF7b': [27756, 27888],
  'ORF8': [27894, 28260],
  'N': [28274, 29534],
  'ORF10': [29558, 29675]
}
# ORF1ab codon position p>4401 <-> ORF1b codon position p-4401

codon = {
  'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
  'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
  'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}
  
# i is nucleobase position from 0
# Return gene name, codon number (from 0), nucleobase pos (from 0) at start of amino acid
def getaa(i):
  i+=1
  for g in genes:
    [a,b]=genes[g]
    if i>=a and i<=b: break
  else: return "????",1,1
  return g,(i-a)//3,i-(i-a)%3-1

# As getaa() but in orf1ab format
def getaa2(i):
  g,c,p=getaa(i)
  if g=='ORF1a': g='orf1ab'
  elif g=='ORF1b': g='orf1ab';c+=4401
  return g,c,p

