import csv,time,calendar,os,json,sys,datetime

# Returns map: {headings} -> [list of entries]
def loadcsv(fn,sep=","):
  with open(fn,'r') as fp:
    reader=csv.reader(fp,delimiter=sep)
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

def datetoday(x):
  if '/' in x:
    t=time.strptime(x+'UTC',"%d/%m/%Y%Z")
  else:
    t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

