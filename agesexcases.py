import time,calendar,os,json,sys,datetime,requests,sys
import numpy as np
from stuff import *

# Publish date or specimen date
usepublishdate=False
discardspecdays=1

# Go back this number of days
ndays=60

# Cache directory
adir='apidata_agesex'

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(3):
    response = requests.get(url+req, timeout=10)
    if response.ok: break
  if not response.ok:
    raise RuntimeError(f'Request failed: { response.text }')
  return response.json()['body'][::-1]


# Convert (eg) string ages '15_19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  if x[-1]=='+': return (int(x[:-1]),150)
  x=x.replace('_to_','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

if usepublishdate:
  print("Using publication date")
  print()
  maxday=apiday()
  # Going up to publish day maxday, which will provide day maxday-1
  minday=maxday-ndays+1
  #print("Date range:",daytodate(minday),"-",daytodate(maxday))
  
  os.makedirs(adir,exist_ok=True)
  l=os.listdir(adir)
  dd=[]
  for day in range(minday-1,maxday+1):
    date=daytodate(day)
    fn=os.path.join(adir,date)
    if os.path.isfile(fn):
      with open(fn,'r') as fp: td=json.load(fp)
    else:
      date=daytodate(day)
      male=get_data('areaType=nation&areaName=England&metric=maleCases&release='+date)
      female=get_data('areaType=nation&areaName=England&metric=femaleCases&release='+date)
      td={**(male[-1]),**(female[-1])}
      del td['metric']
      with open(fn,'w') as fp: json.dump(td,fp,indent=2)
      print("Retrieved api data at",date)
    dd.append(td)
  
  # Put into sensible format
  apiages=[x['age'] for x in dd[-1]['maleCases']]
  ages=sorted([parseage(x) for x in apiages])
  nages=len(ages)
  assert len(dd)==ndays+1
  newcases=np.zeros([ndays,nages,2],dtype=int)
  cumcases=[]
  for day in range(minday-1,maxday+1):
    cumcases.append({a:[0,0] for a in ages})
    for (s,sex) in enumerate(['male','female']):
      for x in dd[day-(minday-1)][sex+'Cases']:
        a=parseage(x['age'])
        v=x['value']
        cumcases[-1][a][s]=v
        if day>=minday: newcases[day-minday,ages.index(a),s]=cumcases[-1][a][s]-cumcases[-2][a][s]
else:
  # Retrieve case data by specimen date
  print("Using specimen date with %d day%s discarded"%(discardspecdays,'' if discardspecdays==1 else 's'))
  print()
  if len(sys.argv)>1:
    release='&release='+sys.argv[1]
  else:
    release=''
  male=get_data('areaType=nation&areaName=England&metric=maleCases'+release)
  female=get_data('areaType=nation&areaName=England&metric=femaleCases'+release)
  male=male[:len(male)-discardspecdays]
  female=female[:len(female)-discardspecdays]
  apiages=[x['age'] for x in male[-1]['maleCases']]
  ages=sorted([parseage(x) for x in apiages])
  nages=len(ages)
  cumcases=np.zeros([ndays+1,nages,2],dtype=int)
  maxday=datetoday(male[-1]['date'])
  minday=maxday-ndays+1
  for (s,sex) in enumerate([male,female]):
    for d in range(ndays+1):
      c=sex[d-ndays-1]
      for x in c[c['metric']]:
        a=ages.index(parseage(x['age']))
        cumcases[d,a,s]=x['value']
  newcases=cumcases[1:]-cumcases[:-1]

gap=7;ave=3
for (s,sex) in enumerate(['male','female','both','male-female','1000*(male-female)/(male+female)']):
  if s<2: cases=newcases[:,:,s]
  elif s==2: cases=newcases.sum(axis=2)
  elif s==3: cases=newcases[:,:,0]-newcases[:,:,1]
  else:
    cases=(newcases[:,:,0]-newcases[:,:,1])/newcases.sum(axis=2)*1000
    totage=newcases.sum(axis=1)
    tot=(totage[:,0]-totage[:,1])/totage.sum(axis=1)*1000
  print('Showing:',sex)
  for ph in range(2-(s>=3)):
    print("     Age:",end='')
    for a in ages: print("%3d  "%a[0],end=' ')
    print("     Total")
    for d in range(ndays):
      print(daytodate(minday+d),end=' ')
      t=cases[d].sum()
      if ph==0: q=1
      else: q=t/1000
      for a in range(nages):
        print("%5d"%((cases[d][a])/q+.5),end=' ')
      if s<4:
        print("= %6d"%(t/q+.5))
      else:
        print("| %6d"%(tot[d]))
    print()
  if s<3:
    print("     Age:",end='')
    for a in ages: print("%3d  "%a[0],end=' ')
    print("     Total")
    for d in range(gap+ave-1,ndays):
      print(daytodate(minday+d),end=' ')
      for a in range(nages):
        print(" %5.2f"%(sum(cases[d-ave+1:d+1,a])/sum(cases[d-gap-ave+1:d-gap+1,a])),end=' ')
      tt=cases.sum(axis=1)
      print("  : %5.2f"%(sum(tt[d-ave+1:d+1])/sum(tt[d-gap-ave+1:d-gap+1])))
    print()


gap=7;ave=3
a0=4;a1=7
ncr=newcases[:,a0:a1,:].sum(axis=1)
cases=(ncr[:,0]-ncr[:,1])/(ncr[:,0]+ncr[:,1])*1000
print('Showing:',sex)
print("     Ages: %d-%d"%(ages[a0][0],ages[a1-1][1]))
for d in range(ndays):
  print(daytodate(minday+d),"%5d"%(cases[d]+.5))

