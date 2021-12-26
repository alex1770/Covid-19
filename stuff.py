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

def datetoday(x):
  format=None
  if '/' in x:
    format='%d/%m/%Y';# E.g., 05/06/2021
  elif ' ' in x:
    format='%d %B %Y';# E.g., 05 June 2021
  elif '-' in x:
    y=x.split('-')
    if len(y[0])==4 and y[1].isdigit():
      format='%Y-%m-%d';# E.g., 2021-06-05
    elif y[1].isalpha() and len(y[2])==2:
      format='%d-%b-%y';# E.g., 29-Oct-20
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

def apiday():
  import pytz
  now=datetime.datetime.now(tz=pytz.timezone('Europe/London'))
  return datetoday(now.strftime('%Y-%m-%d'))-(now.hour<16)

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
  return [s*x for x in adjusted]
    
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
class Date(int):
  def __new__(cls, daydate):
    if type(daydate)==str: return super(cls,cls).__new__(cls,datetoday(daydate))
    else: return super(cls,cls).__new__(cls,daydate)
  def __eq__(self,other): return int(self)==int(Date(other))
  def __hash__(self): return hash(int(Date(self)))
  def __ge__(self,other): return int(self)>=int(Date(other))
  def __le__(self,other): return int(self)<=int(Date(other))
  def __gt__(self,other): return int(self)>int(Date(other))
  def __lt__(self,other): return int(self)<int(Date(other))
  # Date + int = Date
  # Date + str = str
  def __add__(self, other):
    if type(other)==str: return daytodate(int(self))+other
    res = super(Date, self).__add__(other)
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

def getapidata(req):
  import requests
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(5):
    try:
      response = requests.get(url+req, timeout=5)
      if response.ok: break
      error=response.text
    except BaseException as err:
      error=str(err)
  else: raise RuntimeError('Request failed: '+error)
  return response.json()['body'][::-1]

# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

def getcasesbyage(specday,location):
  origages=[(a,a+5) for a in range(0,90,5)]+[(90,150)]
  cachedir=os.path.join(gettopdir(),'apidata_allcaseages')
  if location=='England':
    areatype='nation'
  else:
    areatype='region'
    cachedir+='_'+location
  date=str(Date(specday))
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
    print("Retrieved api data at",date,"for",location)
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
  maxday=Date(maxday)
  astrings=["%d_%d"%a for a in origages]
  reduceages={}
  for (a,astr) in enumerate(astrings):
    for (da,dat) in enumerate(ages):
      if origages[a][0]>=dat[0] and origages[a][1]<=dat[1]: reduceages[astr]=da;break
  nages=len(ages)
  
  # Target save format is
  # filename=publishdate, td[sex][specimendate][agerange] = cumulative cases,
  # having converted agerange to open-closed format and eliminated superfluous ranges, but kept as a string because json can't handle tuples
  # Note that specimendate goes back to the dawn of time, whatever minday is, because we want to save everything.
  # Collect dd[publishdate]=td, td:sex -> specdate -> agestring -> number_of_cases
  
  # Convert to numpy array
  # cc[pastpublishday - (minday-1)][specimenday - (minday-1)][sex 0=m, 1=f][index into ages] = pastpublishday's version of cumulative cases up to specimen day
  # for minday-1 <= pastpublishday <= maxday, minday-1 <= specimenday <= maxday-1
  npub=maxday-minday+1
  nspec=maxday-minday
  cc=np.zeros([npub+1,nspec+1,2,nages],dtype=int)
  smindate=daytodate(minday-1)# Prepare this to compare strings because datetoday is slow
  for pubdate in dd:
    pday=int(pubdate)-(minday-1)
    assert pday>=0
    for sex in dd[pubdate]:
      s=['male','female'].index(sex)
      for specdate in dd[pubdate][sex]:
        if specdate>=smindate:
          sday=datetoday(specdate)-(minday-1)
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

# Return incomplete sample-day correction factors (between 0 and 1) in an array whose
# sample days correspond to (publishday-n, publishday-(n-1), ..., publishday-1).
def getextrap(publishday,location='England'):
  infinity=7
  import numpy as np
  import os,json,sys,datetime

  minday=Date('2021-08-20')
  
  ages=[(a,a+10) for a in range(0,70,10)]+[(70,150)]
  nages=len(ages)
  publishday=Date(publishday)
  
  # Target save format is
  # filename=publishdate, td[sex][specimendate][agerange] = cumulative cases,
  # having converted agerange to open-closed format and eliminated superfluous ranges, but kept as a string because json can't handle tuples
  # Note that specimendate goes back to the dawn of time, whatever minday is, because we want to save everything.
  # Collect dd[publishdate]=td, td:sex -> specdate -> agestring -> number_of_cases
  dd={}
  for day in Daterange(publishday-max(7,infinity-2),publishday+1):
    dd[day]=getcasesbyage(day,location)

  npub,nspec,cc,cn,nn=convcasesbyagetonumpy(dd,minday,publishday,ages=ages)
  gg=cn.sum(axis=2)# Sum over sexes
  
  # Try to undo the effect of delay from specimen to published test result by assuming the pattern is the same as last week's
  # sp[specimenday-minday][age index] = Est no. of samples
  sp=np.zeros([nspec,nages],dtype=float)
  for i in range(nspec):
    if npub-(i+1)>=infinity:
      sp[i]=gg[npub,i,:]
    else:
      #sp[i]=gg[npub,i,:]/gg[npub-7,i-7,:]*gg[i-7+infinity+1,i-7,:]
      base=np.array([gg[npub-r,i-r] for r in range(max(8,infinity-(npub-(i+1))))])
      targ7=gg[i+infinity-7,i-7]
      f0=1+sum((gg[npub,i-r]-gg[npub-1,i-r])/base[r] for r in range(1,infinity-(npub-(i+1))))
      f1=targ7/base[7]
      sp[i]=base[0]*(0.55*f0+0.45*f1)
    
  return gg[npub,:,:].sum(axis=1)/sp.sum(axis=1)

#def getextrap(publishday,location='England',minday=Date('2021-08-20')):
#  sp=getcorrectedcases(minday,publishday,location=location)
  
