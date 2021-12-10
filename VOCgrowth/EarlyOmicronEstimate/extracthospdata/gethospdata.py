import os,json

# Decompile hospital surveillance reports from NICD, South Africa and put into json format
# From https://www.nicd.ac.za/diseases-a-z-index/disease-index-covid-19/surveillance-reports/daily-hospital-surveillance-datcov-report/

datafile='SouthAfricaHospData.json'
textdir='textpdfs'

pdfs=os.listdir(textdir)
pdfs.sort(reverse=True)

if os.path.isfile(datafile):
  with open(datafile) as fp: data=json.load(fp)
else:
  data={}

provinces=['Eastern Cape','Free State','Gauteng','KwaZulu-Natal','Limpopo','Mpumalanga','North West','Northern Cape','Western Cape']
headings=['Facilities Reporting','Admissions to Date','Died to Date','Discharged to Date','Currently Admitted','Currently in ICU','Currently Ventilated','Currently Oxygenated','Admissions in Previous Day']
keyd={x.split()[0]:(len(x.split()),x) for x in provinces}
keyd['Total']=(1,'South Africa')
for date in pdfs:
  with open(os.path.join(textdir,date)) as fp: text=fp.read()
  text2=text.split('\n')

  # Main table (National version)
  i0=1000000
  intable=False
  tabletitle='Summary of reported COVID-19 admissions by province'
  if date not in data: data[date]={}
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

  # Admissions & deaths by age group
  # Before 2021-10-27 there was a different format with missing elements that I'm not going to try to parse
  #if date>='2021-10-27':
    

# Order by reverse date to convenience of (human) reading file
odata={x:data[x] for x in sorted(list(data),reverse=True)}

with open(datafile,'w') as fp: json.dump(odata,fp,indent=2)
