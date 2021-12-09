import os,bs4,datetime,subprocess,json

datafile='SouthAfricaHospData.json'

# pdf filenames don't always contain their date, so get these from the index page
path2date={}
with open('www.nicd.ac.za/index.html?p=25308') as fp:
  soup=bs4.BeautifulSoup(fp,'html5lib')
for x in soup.find_all('a'):
  href=x.attrs['href']
  text=x.text
  if href[-4:]=='.pdf' and 'NICD COVID-19 SURVEILLANCE IN SELECTED HOSPITALS' in text:
    f0=text.find('(');f1=text.find(')')
    t1=text[f0+1:f1].split()
    t1[1]=t1[1][:3]
    date=datetime.datetime.strptime(' '.join(t1),'%d %b %Y').strftime('%Y-%m-%d')
    path='www.nicd.ac.za'+href[22:]
    path2date[path]=date

pdfs=[]
for x in os.walk('www.nicd.ac.za/wp-content/uploads'):
  for y in x[2]:
    if y[-4:]=='.pdf' and 'Weekly' not in y and 'Monthly' not in y:
      path=x[0]+'/'+y
      if path in path2date: date=path2date[path]
      else: date=path[-12:-8]+'-'+path[-8:-6]+'-'+path[-6:-4]
      pdfs.append((date,path))
pdfs.sort(reverse=True)


data={}
provinces=['Eastern Cape','Free State','Gauteng','KwaZulu-Natal','Limpopo','Mpumalanga','North West','Northern Cape','Western Cape']
headings=['Facilities Reporting','Admissions to Date','Died to Date','Discharged to Date','Currently Admitted','Currently in ICU','Currently Ventilated','Currently Oxygenated','Admissions in Previous Day']
keyd={x.split()[0]:(len(x.split()),x) for x in provinces}
keyd['Total']=(1,'South Africa')
for (date,path) in pdfs:
  # pdftotext -layout <file> -
  po=subprocess.Popen(['/usr/bin/pdftotext','-layout',path,'-'],stdout=subprocess.PIPE,encoding='utf-8')
  p=po.stdout
  text=p.read()
  p.close()
  po.wait()
  
  text2=text.split('\n')
  i0=1000000
  print(date)
  intable=False
  tabletitle='Summary of reported COVID-19 admissions by province'
  data[date]={}
  for (i,l) in enumerate(text2):
    l=l.strip()
    if intable:
      ll=l.split()
      if len(ll)==0: continue
      if ll[0] in keyd:
        (n,location)=keyd[ll[0]]
        data[date][location]={}
        for desc,num in zip(headings,ll[n:]):data[date][location][desc]=int(num)
      if ll[0]=='Total': break
    else:
      if l[:len(tabletitle)]==tabletitle and 'sector' in l: intable=True

with open(datafile,'w') as fp: json.dump(data,fp,indent=2)
