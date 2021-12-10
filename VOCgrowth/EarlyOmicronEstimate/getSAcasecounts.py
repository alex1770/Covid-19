from stuff import *
import requests,bs4,pytz,datetime,csv,sys

mindate='2021-06-17'
if len(sys.argv)>1: mindate=sys.argv[1]
datafile='SAcasecounts.csv'

provinces=['Total','Eastern Cape', 'Free State', 'Gauteng', 'KwaZulu-Natal', 'Limpopo', 'Mpumalanga', 'North West', 'Northern Cape', 'Western Cape']

d=datetime.datetime.now(pytz.timezone("Africa/Johannesburg"))
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<12: today-=1# Don't look for new cases before noon local time

newdict={}    
headings=['Date']+provinces
if os.path.isfile(datafile):
  with open(datafile,'r') as fp:
    r=csv.reader(fp)
    h=next(r)
    if h==headings:
      for crow in r:
        newdict[datetoday(crow[0])]=[int(x) for x in crow[1:]]

totdict={}# Not currently used
minday=datetoday(mindate)
for day in range(minday,today+1):
  date=daytodate(day)
  if day in newdict: continue
  d=datetime.datetime.strptime(date,'%Y-%m-%d')
  if date=='2021-10-30':
    # For some reason, the normal URL doesn't exist for 30 October 2021
    url='https://www.nicd.ac.za/36074-2/'
  elif (date>='2021-07-01' and date<='2021-07-06') or date=='2021-07-09' or (date>='2021-09-06' and date<='2021-09-09'):
    # It uses different day convention for these days
    suffix=d.strftime('%d-%B-%Y').strip().lower()
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-'+suffix
  else:
    suffix=d.strftime('%e-%B-%Y').strip().lower()
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-'+suffix
  print('Loading',date)
  resp=requests.get(url)
  if not resp.ok:
    #if day<today: raise RuntimeError('Could not load '+url)
    #break
    if day<today: print('Could not load '+url,file=sys.stderr)
    continue
  f=max(resp.text.find('PROVINCIAL BREAKDOWN'),0)
  soup=bs4.BeautifulSoup(resp.text[f:],'html5lib')
  tab=soup.find('table')
  ok=1
  col0=col1=col2=None# Before a certain date, there wasn't a 'new cases' column
  newcases={}
  totcases={}
  rownum=0
  for row in tab.find_all('tr'):
    if not ok: break
    if rownum>0 and col0==None:
      # Some tables, e.g., 2021-06-02, don't have a province column, so assume default order
      if rownum<len(provinces): prov=provinces[rownum]
      else: prov='Total'# It's possible this will pick up the 'Unknown' row, but it will later be overwritten by the correct 'Total' row
    colnum=0
    for col in row.find_all('td'):
      if 'width' in col.attrs and int(col.attrs['width'])<20: continue
      t=col.text.strip()
      if t=='': ok=0;break
      if rownum==0:
        if t=='Province': col0=colnum
        if 'new cases on' in t.lower(): col1=colnum
        if 'total cases for' in t.lower(): col2=colnum
      else:
        if colnum==col0: prov=t.strip('*')
        if colnum==col1: newcases[prov]=int(t.replace(',','').replace(' ',''))
        if colnum==col2: totcases[prov]=int(t.replace(',','').replace(' ',''))
      colnum+=1
    rownum+=1
  if col1!=None: newdict[day]=[newcases[prov] for prov in provinces]
  totdict[day]=[totcases[prov] for prov in provinces]

with open(datafile,'w') as fp:
  w=csv.writer(fp)
  w.writerow(headings)
  for day in sorted(list(newdict)):
    w.writerow([daytodate(day)]+newdict[day])
