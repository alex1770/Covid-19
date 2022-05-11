from stuff import *
import requests,bs4,pytz,datetime,csv,sys

mindate='2020-12-12'
if len(sys.argv)>1: mindate=sys.argv[1]
datafile='SAcasecounts.csv'

provinces=['Eastern Cape', 'Free State', 'Gauteng', 'KwaZulu-Natal', 'Limpopo', 'Mpumalanga', 'North West', 'Northern Cape', 'Western Cape']

d=datetime.datetime.now(pytz.timezone("Africa/Johannesburg"))
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<12: today-=1# Don't look for new cases before noon local time

newdict={}    
headings=['Date','Total']+provinces
if os.path.isfile(datafile):
  with open(datafile,'r') as fp:
    r=csv.reader(fp)
    h=next(r)
    if h==headings:
      for crow in r:
        newdict[datetoday(crow[0])]=[int(x) for x in crow[2:]]

def getint(t):
  if '\xa0\xa0' in t: t=t[:t.find('\xa0\xa0')]# Sometimes the next column is concatenated in the same cell, separated by non-breaking spaces
  return int(t.replace(',','').replace(' ','').replace('\xa0',''))
      
totdict={}
minday=datetoday(mindate)
for day in range(minday,today+1):
  date=daytodate(day)
  if day in newdict and day+1 in newdict: continue
  d=datetime.datetime.strptime(date,'%Y-%m-%d')
  # Special case URLs
  if date=='2020-12-31':
    #url='https://sacoronavirus.co.za/2020/12/31/update-on-covid-19-31st-december-2020/'
    totdict[day]=[170687, 62670, 287018, 199983, 25076, 36707, 40185, 25314, 209521]
    continue
  elif date=='2021-01-06':
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-06-january-20210/'
  elif date=='2021-01-20':
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-19-jan-2021-2/'
  elif date=='2021-05-25':
    # This date is missing - can't find an alternative
    continue
  elif date=='2021-06-08':
    # 8 June 2021 doesn't exist in the records, but can infer from 9 June
    totdict[day]=[1712939-8881, 199165-253, 105411-620, 484116-5111, 341871-537, 66642-229, 84510-363, 79982-650, 53433-335, 297809-783]
    continue
  elif date=='2021-05-08':
    url='https://www.nicd.ac.za/28717-2/'
  elif date=='2021-06-11':
    url='https://www.nicd.ac.za/29584-2/'
  elif date=='2021-10-30':
    url='https://www.nicd.ac.za/36074-2/'
  elif date=='2021-12-19':
    url='https://www.nicd.ac.za/37362-2/'
  elif date=='2021-12-29':
    url='https://www.nicd.ac.za/confirmed-cases-of-covid-19-in-south-africa-29-december-2021/'
  elif date=='2021-12-30':
    url='https://www.nicd.ac.za/confirmed-cases-of-covid-19-in-south-africa-30-december-2021/'
  elif date<'2021-01-10':
    suffix=d.strftime('%d-%b-%Y').strip().lower()
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-'+suffix
  elif (date>='2021-07-01' and date<='2021-07-06') or date=='2021-07-09' or (date>='2021-09-06' and date<='2021-09-09') or (date>='2022-01-06' and date<='2022-01-09'):
    # E.g. https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-06-january-2022
    suffix=d.strftime('%d-%B-%Y').strip().lower()
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-'+suffix
  elif date<'2021-04-01' or (date>='2021-04-10' and date<'2021-05-01'):
    suffix=d.strftime('%e-%b-%Y').strip().lower()
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-'+suffix
  elif date=='2022-01-17':
    url='https://www.nicd.ac.za/confirmed-cases-of-covid-19-in-south-africa-17-january-2022/'
  elif date=='2022-02-07':
    url='https://www.nicd.ac.za/confirmed-cases-of-covid-19-in-south-africa-7-february-2022/'
  else:
    suffix=d.strftime('%e-%B-%Y').strip().lower()
    url='https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-'+suffix
  print('Loading',date,'from',url)
  resp=requests.get(url)
  if not resp.ok:
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
      if rownum-1<len(provinces): prov=provinces[rownum-1]
      else: break
    colnum=0
    for col in row.find_all('td'):
      if 'width' in col.attrs and int(col.attrs['width'])<20: continue
      t=col.text.strip()
      if t=='': break#ok=0;break# Can't do ok=0 for 2021-06-16. Can't remember why needed ok=0
      if rownum==0:
        if t=='Province': col0=colnum
        if 'new cases on' in t.lower(): col1=colnum
        if 'total cases for' in t.lower() or (date<'2021-05-01' and t[:6]=='Cases '): col2=colnum
      else:
        if colnum==col0: prov=t.strip('*')
        if colnum==col1: newcases[prov]=getint(t)
        if colnum==col2: totcases[prov]=getint(t)
      colnum+=1
    rownum+=1
  if col1!=None: newdict[day]=[newcases[prov] for prov in provinces]
  totdict[day]=[totcases[prov] for prov in provinces]

# Interpolate missing 2021-05-25 data
# https://www.nicd.ac.za/latest-confirmed-cases-of-covid-19-in-south-africa-25-may-2021/ only gives 24 May data + overall total
if 18771 in totdict and 18773 in totdict and 18772 not in totdict:
  tot18771=sum(totdict[18771])
  tot18772=1640932# Total reported on 25 May
  tot18773=sum(totdict[18773])
  al=(tot18772-tot18771)/(tot18773-tot18771)
  totdict[18772]=[int(x+al*(y-x)+.45) for (x,y) in zip(totdict[18771],totdict[18773])]
    
for day in range(minday+1,today+1):
  if day not in newdict and day in totdict and day-1 in totdict:
    newdict[day]=[d1-d0 for (d0,d1) in zip(totdict[day-1],totdict[day])]

with open(datafile,'w') as fp:
  w=csv.writer(fp)
  w.writerow(headings)
  for day in sorted(list(newdict)):
    w.writerow([daytodate(day),sum(newdict[day])]+newdict[day])
