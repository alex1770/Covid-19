import sys,re,bs4,time

website='https://www.worldometers.info/coronavirus/'
from urllib.request import Request, urlopen
req = Request(website, headers={'User-Agent': 'Mozilla/5.0'})
u = urlopen(req)
soup=bs4.BeautifulSoup(u,"html5lib")
mintime="16-00-00"

# 2020-02-29 -> 2020-03-01 etc
def tomorrow(x):
  (y,m,d)=map(int,x.split('-'))
  d+=1
  if d>[31,28,31,30,31,30,31,31,30,31,30,31][m-1]+(m==2 and y%4==0):
    d=1;m+=1
    if m==13:
      m=1;y+=1
  return "%04d-%02d-%02d"%(y,m,d)
    
# Convert e.g., "Feb 14" to 2020-02-14
def isofy(x):
  y=x.split()
  if len(y)!=2: return "<unrecogniseddate>"
  f="JanFebMarAprMayJunJulAugSepOctNovDec".find(y[0])
  if f<0 or not y[1].isdigit(): return "<unrecogniseddate>"
  return "2020-%02d-%02d"%(f//3+1,int(y[1]))

keys=["'coronavirus-cases-linear'",
      "'coronavirus-deaths-linear'",
      "<recovered>",
      "'graph-active-cases-total'"
      ]

def pr(x):
  if x>=0: return str(x)
  return '?'

def prsub(x,y):
  if x>=0 and y>=0 and x>=y: return str(x-y)
  return '?'

# Ignore countries with less than 'mincount' confirmed cases
mincount=100

print("date,country,new_cases,new_deaths,total_cases,total_deaths,recovered,active")
tab=soup.find(id="main_table_countries_yesterday")
lastcount=mincount
for row in tab.find_all("tr"):
  for col in row.find_all('td'):
    if col!=None:
      x=col.find("a")
      if x!=None:
        lastcount=0
        country=x.text
        c=x.attrs['href']
        if c[:8]=='country/':
          #print("%-15s"%country,website+c,file=sys.stderr)
          print("Parsing",country,file=sys.stderr)
          req2=Request(website+c, headers={'User-Agent': 'Mozilla/5.0'})
          u2=urlopen(req2)
          soup2=bs4.BeautifulSoup(u2,"html5lib")
          d={}
          for x in soup2.find_all('script'):
            text=x.text.replace('\n',' ')
            r=re.search(r"Highcharts.chart\((.*?),",text)
            if r!=None:
              name=r.group(1)
              r=re.search(r"xAxis:\s*{\s*categories:\s*\[(.*?)\s*\]\s*}", text)
              if r!=None:
                dates=[isofy(x.strip('"')) for x in r.group(1).strip().split(",")]
              else: dates=[]
              r=re.search(r"series.*data:\s*\[(.*?)\]\s*}", text)
              if r!=None: values=r.group(1).strip().split(",")
              else: values=[]
              if name in keys:
                p=keys.index(name)
                for (dt,vl) in zip(dates,values):
                  if dt not in d: d[dt]=[-1]*4
                  d[dt][p]=int(vl)
              #assert len(dates)==len(values)
          for dt in d:
            if d[dt][0]>=0 and d[dt][1]>=0 and d[dt][3]>=0: d[dt][2]=d[dt][0]-d[dt][1]-d[dt][3]
          # Possibly add one more entry from the "Latest Updates" section, but be cautious
          # because sometimes the entries are premature and too low.
          if len(dates)>=2:
            prev=dates[-2]
            last=dates[-1]
            targ=tomorrow(last)
            now=time.strftime("%Y-%m-%d-%H-%M-%S",time.gmtime())
            if now>=targ+"-"+mintime:
              newc=newd=-1
              for x in soup2.find_all('div'):
                if 'id' in x.attrs and x.attrs["id"]=="newsdate"+targ:
                  #print(x.text)
                  y=x.text.split()
                  if "deaths" in y:
                    i=y.index("deaths")
                    if i>=2 and y[i-1]=="new":
                      r=y[i-2].replace(',','')
                      if r.isdigit(): newd=int(r)
                  if "cases" in y:
                    i=y.index("cases")
                    if i>=2 and y[i-1]=="new":
                      r=y[i-2].replace(',','')
                      if r.isdigit(): newc=int(r)
              # Reject this "Latest Updates" entry if not plausibly big enough:
              oldc=d[last][0]-d[prev][0]
              oldd=d[last][1]-d[prev][1]
              if (newc>=0 and newc>=(d[last][0]-d[prev][0])/2) and (newd>=0 and newd>=(d[last][1]-d[prev][1])/2):
                d[targ]=[d[last][0]+newc,d[last][1]+newd,-1,-1]
          if d!={}:
            lastcount=d[sorted(list(d))[-1]][0]
            #if lastcount<mincount: break
          prc=prd=0
          for dt in sorted(list(d)):
            print("%s,\"%s\",%s,%s,%s,%s,%s,%s"%(dt,country,prsub(d[dt][0],prc),prsub(d[dt][1],prd),pr(d[dt][0]),pr(d[dt][1]),pr(d[dt][2]),pr(d[dt][3])))
            prc,prd=d[dt][:2]
      #if lastcount<mincount: break
