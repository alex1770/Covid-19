import csv,time,calendar,os,json,sys,datetime,requests,pytz

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
      out[h]=[int(x) for x in out[h]]
    except ValueError:
      try:
        out[h]=[float(x) for x in out[h]]
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
