from stuff import *
import json,requests,csv,os
import numpy as np
from math import sqrt

np.set_printoptions(precision=4,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=10000)

# ONS 2020 pop ages originally given as 0-18, 18-25, 25-30, ..., 75-80, 80+. Interpolate children age bands
# Get NIMSpop from dashboard (vaccinationsAgeDemographics metric)
ONSages=[(0,12),         (12,16),       (16,18),       (18,25), (25,30), (30,35), (35,40), (40,45), (45,50), (50,55), (55,60), (60,65), (65,70), (70,75), (75,80), (80,150)]
ONSpop= [12093288*12/18, 12093288*4/18, 12093288*2/18, 4709589, 3771493, 3824652, 3738209, 3476303, 3638639, 3875351, 3761782, 3196813, 2784300, 2814128, 2009992, 2855599]

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(10):
    try:
      response = requests.get(url+req, timeout=5)
      if response.ok: break
      error=response.text
    except BaseException as err:
      error=str(err)
  else: raise RuntimeError('Request failed: '+error)
  return response.json()['body'][::-1]

# Convert (eg) string ages '15-19', '15_to_19', '60+' to (15,20), (15,20), (60,150) respectively
def parseage(x):
  x=x.strip()
  if x[-1]=='+': return (int(x[:-1]),150)
  if x[:6]=='Under ': return (0,int(x[6:]))
  x=x.replace('_to_','_').replace('-','_')# cater for 65_to_69 and 65_69 formats
  aa=[int(y) for y in x.split("_")]
  return (aa[0],aa[1]+1)

# Output in muggle format
def unparseage(ages):
  (x,y)=ages
  if y<150: return "%d-%d"%(x,y-1)
  return "%d+"%x

rawvax=get_data('areaType=nation&areaName=England&metric=vaccinationsAgeDemographics')

#pmd=loadcsv("vaccine-surveillance-reports/vaccine surveillance.csv")
#ukhsaages=sorted(parseage(a) for a in set(pmd['Age']))
#ag=(40,50)
#ty='Cases'
#for (a,d,t,ul,d0,dhalf,d1,d2) in zip(pmd['Age'],pmd['Date'],pmd['Type'],pmd['Unlinked'],
#                                     pmd['Not vaccinated'],
#                                     pmd['Received one dose 1-20 days before specimen date'],
#                                     pmd['Received one dose, 21 days or more before specimen date'],
#                                     pmd['Second dose 14 days or more before specimen date']):
#  if t!=ty: continue
#  if parseage(a)!=ag: continue
#  #print(a,d,t,ul,d0,dhalf,d1,d2)


# Table 1g is from ONS infection survey data at https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk
# It's a csv page saved from a spreadsheet with messy formatting, so have to parse manually here.
# Make survey[][][][]:
# survey[ageindex][doses-1][dateindex]=[ONS modelled prob, ONS prob low, ONS prob high, number reported receiving >=doses, number in sample]
days=[]# Week commencing
with open('ONSinfectionsurvey.table1g.csv','r') as fp:
  reader=csv.reader(fp)
  surveyages=[];agecols=[]
  while 1:
    row=next(reader)
    for (i,a) in enumerate(row):
      if a[:4]=='Age ':
        surveyages.append(parseage(a[4:]))
        agecols.append(i)
    if surveyages: break
  nsurveyages=len(surveyages)
  survey=[[[],[]] for i in range(nsurveyages)]
  headings=next(reader)
  parseheadings={}
  for (i,c) in enumerate(agecols):
    for j in range(c,c+agecols[1]-agecols[0]):
      h=headings[j]
      if '1 or more' in h: d=0
      elif '2 or more' in h: d=1
      if 'Modelled' in h: parseheadings[(i,d,0)]=j
      elif 'Lower' in h: parseheadings[(i,d,1)]=j
      elif 'Upper' in h: parseheadings[(i,d,2)]=j
      elif 'Number' in h and 'receiv' in h: parseheadings[(i,d,3)]=j
      elif 'adults in sample' in h: parseheadings[(i,d,4)]=j
  while 1:
    row=next(reader)
    if not(row[0][:1].isdigit() and ' to ' in row[0]): break
    dd=[datetoday(time.strftime("%Y-%m-%d",time.strptime(x,"%d %B %Y"))) for x in row[0].split(' to ')]
    days.append(dd[0])
    assert dd[1]-dd[0]==6# If this changes then would need to change averaging below
    for (a,d,i) in parseheadings:
      c=parseheadings[a,d,i]
      if i==0: survey[a][d].append([-1]*5)
      if row[c]!='-': survey[a][d][-1][i]=float(row[c].replace(',',''))/(100 if i<3 else 1)
survey=np.array(survey)

ndates=len(days)

offset=0

# Now parse rawvax to make:
# vax[ageindex][doses-1][dateindex] = number vaccinated
#
datetoindex={}
for (i,day) in enumerate(days):
  for d in range(day,day+7):
    datetoindex[daytodate(d+offset)]=i

numvaxages=sorted([parseage(v['age']) for v in rawvax[0][rawvax[0]['metric']]])
numvaxages2surveyages={}
for (a,b) in numvaxages:
  for (i,(c,d)) in enumerate(surveyages):
    if max(a,c)<min(b,d):
      assert a>=c and b<=d
      numvaxages2surveyages[(a,b)]=i;break

ONSpop_surveyages=[0]*nsurveyages
for (i,(a,b)) in enumerate(ONSages):
  for (j,(c,d)) in enumerate(surveyages):
    if max(a,c)<min(b,d):
      assert a>=c and b<=d
      ONSpop_surveyages[j]+=ONSpop[i]

# The NIMS populations seem to be given as the same for all dates, so use the last one here
vv=rawvax[-1][rawvax[-1]['metric']]
NIMSpop_surveyages=[0]*nsurveyages
for v in vv:
  (a,b)=parseage(v['age'])
  for (i,(c,d)) in enumerate(surveyages):
    if max(a,c)<min(b,d):
      assert a>=c and b<=d
      NIMSpop_surveyages[i]+=int(v['VaccineRegisterPopulationByVaccinationDate'])


vax=np.zeros([nsurveyages,2,ndates])
for vv in rawvax:
  if vv['date'] not in datetoindex: continue
  d=datetoindex[vv['date']]
  for v in vv[vv['metric']]:
    age=parseage(v['age'])
    if age not in numvaxages2surveyages: continue
    a=numvaxages2surveyages[age]
    vax[a][0][d]+=v['cumPeopleVaccinatedFirstDoseByVaccinationDate']
    vax[a][1][d]+=v['cumPeopleVaccinatedSecondDoseByVaccinationDate']
vax/=7

graphdir='graphs'
os.makedirs(graphdir,exist_ok=True)
for a in range(nsurveyages):
  print()
  print("Age range:",unparseage(surveyages[a]))
  print("All populations are in millions")
  print("Week/comm      ONS est, dose >=1          ONS est, dose >=2           ONS pop  NIMS pop")
  print("=========      ======================     ======================      =======  ========")
  for d in range(ndates):
    print(daytodate(days[d]),end='')
    for doseind in range(2):
      su=survey[a][doseind]
      vv=vax[a][doseind]
      if su[d][0]>0: print("    %5.2f ( %5.2f - %5.2f )"%(vv[d]/su[d][0]/1e6, vv[d]/su[d][2]/1e6, vv[d]/su[d][1]/1e6),end='')
      else: print("        - (     - -     - )",end='')
    print("      %7.2f"%(ONSpop_surveyages[a]/1e6),end='')
    print("   %7.2f"%(NIMSpop_surveyages[a]/1e6),end='')
    print()
  print()


  
  data=[]
  maxpop=0
  for doseind in range(2):
    col=['"green"','"blue"'][doseind]
    ep=vax[a][doseind][:,None]/survey[a][doseind]/1e6# Estimated population
    title='(Number vaccinated with ≥%d dose%s) ÷ (Survey estimate of proportion so-vaccinated).'%(doseind+1,doseind*'s')
    e=ep[-1]
    title+='   Last datapoint: %.2fm (%.2fm - %.2fm)'%(e[0],e[1],e[2])
    data.append({
      'with': ('filledcurves',2),
      'title': '',
      'values': [(daytodate(day+3),e[2],e[1]) for (day,e) in zip(days,ep) if e[0]>0],
      'extra': 'lc '+col
    })
    maxpop=max(maxpop,max(e[1] for e in ep))
    data.append({
      'title': title,
      'values': [(daytodate(day+3),e[0]) for (day,e) in zip(days,ep) if e[0]>0],
      'extra': 'lc '+col
    })
  ons=ONSpop_surveyages[a]/1e6
  data.append({
    'title': 'ONS 2020 population (%.2fm)'%ons,
    'values': [(daytodate(days[0]+3),ons),(daytodate(days[-1]+3),ons)],
    'extra': 'lc 4'
  })
  nims=NIMSpop_surveyages[a]/1e6
  data.append({
    'title': 'NIMS population (%.2fm)'%nims,
    'values': [(daytodate(days[0]+3),nims),(daytodate(days[-1]+3),nims)],
    'extra': 'lc 1'
  })
  maxpop=max(maxpop,ons,nims)
  ar=unparseage(surveyages[a])
  title='Estimated population of age band %s using vaccination survey\\nData sources: ONS antibody and vaccination survey, UKHSA dashboard'%ar
  makegraph(title=title, data=data, ylabel='Estimated population (millions)', outfn=os.path.join(graphdir,'PopEst%s.png'%ar),
            extra=['set key top left','set style fill transparent solid 0.25'],ranges='[:] [:%f]'%(maxpop*1.15),interval=86400*14)


  
  data=[]
  for doseind in range(2):
    col=['"green"','"blue"'][doseind]
    vv=vax[a][doseind]/vax[a][doseind][-1]
    su=survey[a][doseind]/survey[a][doseind][-1][0]
    data.append({
      'with': ('filledcurves',2),
      'title': '',
      'values': [(daytodate(day+3),e[1],e[2]) for (day,e) in zip(days,su) if e[0]>0],
      'extra': 'lc '+col
    })
    data.append({
      'title': '(Survey estimate of proportion vaccinated with ≥%d dose%s) ÷ (Latest estimate of same)'%(doseind+1,doseind*'s'),
      'values': [(daytodate(day+3),e[0]) for (day,e) in zip(days,su) if e[0]>0],
      'extra': 'lc '+col
    })
    col=['"dark-olivegreen"','"black"'][doseind]
    data.append({
      'title': '(Number vaccinated with ≥%d dose%s) ÷ (Latest number of same)'%(doseind+1,doseind*'s'),
      'values': [(daytodate(day+3),v) for (day,v) in zip(days,vv)],
      'extra': 'lc '+col
    })
  ar=unparseage(surveyages[a])
  title='Comparison of growth of vaccinated %s year olds: actual numbers vs vaccine survey (artificially normalised to end at 1)\\nData sources: ONS antibody and vaccination survey, UKHSA dashboard'%ar
  makegraph(title=title, data=data, ylabel='Proportion of latest value', outfn=os.path.join(graphdir,'GrowthComparison%s.png'%ar),
            extra=["set key top left",'set style fill transparent solid 0.25'], interval=86400*14)


  
  data=[]
  vv=vax[a][1]/vax[a][0]
  su=survey[a][1][:,0]/survey[a][0][:,0]
  data.append({
    'title': '(Survey estimate of proportion vaccinated with ≥2 doses) ÷ (Survey estimate with ≥1 dose)',
    'values': [(daytodate(day+3),e) for (day,e) in zip(days,su) if e>0],
    'extra': 'lc "blue"'
  })
  data.append({
    'title': '(Number vaccinated with ≥2 doses) ÷ (Number with ≥1 dose)',
    'values': [(daytodate(day+3),v) for (day,v) in zip(days,vv)],
    'extra': 'lc "green"'
  })
  ar=unparseage(surveyages[a])
  title='Proportion of vaccinated %s year olds who have had ≥2 doses, comparing survey estimate with actual numbers\\nData sources: ONS antibody and vaccination survey, UKHSA dashboard'%ar
  makegraph(title=title, data=data, ylabel='Ratio of ≥2 doses to ≥1 dose', outfn=os.path.join(graphdir,'2ndDoseComparison%s.png'%ar),
            extra=["set key top left",'set style fill transparent solid 0.25'],ranges='[:] [0:1.1]', interval=86400*14)
  
