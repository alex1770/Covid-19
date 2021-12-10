import os,bs4,datetime,subprocess

# Convert hospital surveillance reports from NICD, South Africa from pdf to text file, with consistent date naming
# From https://www.nicd.ac.za/diseases-a-z-index/disease-index-covid-19/surveillance-reports/daily-hospital-surveillance-datcov-report/

textdir='textpdfs'

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

os.makedirs(textdir,exist_ok=True)
for (date,path) in pdfs:
  dest=os.path.join(textdir,date)
  if not os.path.isfile(dest):
    print(date)
    # pdftotext -layout <file> -
    po=subprocess.Popen(['/usr/bin/pdftotext','-layout',path,'-'],stdout=subprocess.PIPE,encoding='utf-8')
    p=po.stdout
    text=p.read()
    p.close()
    po.wait()
    with open(dest,'w') as fp:
      fp.write(text)
