from stuff import *
import requests,bs4,pytz,datetime,csv,sys

mindate='2021-11-01'
if len(sys.argv)>1: mindate=sys.argv[1]
datafile='SAcasecounts.csv'

provinces=['Total','Eastern Cape', 'Free State', 'Gauteng', 'KwaZulu-Natal', 'Limpopo', 'Mpumalanga', 'North West', 'Northern Cape', 'Western Cape']

d=datetime.datetime.now(pytz.timezone("Africa/Johannesburg"))
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<12: today-=1# Don't look for new cases before noon

csvdict={}    
headings=['Date']+provinces
if os.path.isfile(datafile):
  with open(datafile,'r') as fp:
    r=csv.reader(fp)
    h=next(r)
    if h==headings:
      for crow in r:
        csvdict[crow[0]]=crow

for day in range(datetoday(mindate),today+1):
  date=daytodate(day)
  if date in csvdict: continue
  d=datetime.datetime.strptime(date,'%Y-%m-%d')
  if date=='2021-10-30':
    # For some reason, the normal URL doesn't exist for 30 October 2021
    url='https://www.nicd.ac.za/36074-2/'
  else:
    suffix=d.strftime('%e-%B-%Y').strip().lower()
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-'+suffix
  print('Loading',date)
  resp=requests.get(url)
  if not resp.ok:
    if day<today: raise RuntimeError('Could not load '+url+suffix)
    break
  f=resp.text.find('PROVINCIAL BREAKDOWN')
  soup=bs4.BeautifulSoup(resp.text[f:],'html5lib')
  tab=soup.find('table')
  first=1
  dat={}
  for row in tab.find_all('tr'):
    i=0
    for col in row.find_all('td'):
      t=col.text.strip()
      if first:
        if t=='Province': col0=i
        if 'new cases on' in t.lower(): col1=i
      else:
        if i==col0: prov=t.strip('*')
        if i==col1: cases=int(t.replace(',','').replace(' ',''))
      i+=1
    if first:
      first=0
    else:
      dat[prov]=cases
  crow=[date]+[dat[prov] for prov in provinces]
  csvdict[date]=crow

with open(datafile,'w') as fp:
  w=csv.writer(fp)
  w.writerow(headings)
  for date in sorted(list(csvdict)):
    w.writerow(csvdict[date])
