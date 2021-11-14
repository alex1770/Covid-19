import csv,time,calendar,os,json,sys,datetime,requests,pytz,subprocess

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
def csvrows_it(it,reqheadings,sep=","):
  reader=csv.reader(it,delimiter=sep)
  headings=next(reader)
  dec={}
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
    yield [row[i] for i in cols]
  
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
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  response=requests.get(url+req,timeout=10)
  if not response.ok:
    raise RuntimeError('Request failed: '+response.text)
  return response

def apiday():
  now=datetime.datetime.now(tz=pytz.timezone('Europe/London'))
  return datetoday(now.strftime('%Y-%m-%d'))-(now.hour<16)

def makegraph(title='A graph', data=[], mindate='0000-00-00', ylabel='', outfn='temp.png', extra=[]):
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
