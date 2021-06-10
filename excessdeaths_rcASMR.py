# Trying to replicate (and possibly improve upon) rcASMR, rASMR calculations from bottom of https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/articles/comparisonsofallcausemortalitybetweeneuropeancountriesandregions/januarytojune2020
# which uses the Institute and Faculty of Actuaries' Continuous Mortality Investigation method as described in working paper 111 here: https://www.actuaries.org.uk/learn-and-develop/continuous-mortality-investigation/cmi-working-papers/mortality-projections/cmi-working-paper-111

# Mortality source: https://data.europa.eu/euodp/en/data/dataset/QrtzdXsI5w26vnr54SIzpQ giving demo_r_mwk_05.tsv
# Explanation: https://ec.europa.eu/eurostat/cache/metadata/en/demomwk_esms.htm
# Ignoring week 99, which appears to contain no data in the examples I've looked at.
# The week number (1-53) is the ISO 8601 week: weeks run Mon-Sun and are assigned to the year in which Thursday falls.

# For population source, can use either Eurostat or UN WPP
#
# Eurostat:
# https://ec.europa.eu/eurostat/web/products-datasets/product?code=urt_pjangrp3 giving urt_pjangrp3.tsv
# (https://ec.europa.eu/eurostat/en/web/products-datasets/-/demo_pjangroup goes back further but doesn't have a 85-89 age group)
# and for population projections, use
# EU: https://ec.europa.eu/eurostat/web/products-datasets/product?code=proj_19np
# UK: (uses midpoints of years): https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/tablea21principalprojectionukpopulationinagegroups
#     or                         https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/tablel21oldagestructurevariantukpopulationinagegroups
# Not currently using population projections for UK (so UK 2020, 2021 projections are done by simple linear interpolation)
#
# UN WPP:
# https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2019_POP_F15_1_ANNUAL_POPULATION_BY_AGE_BOTH_SEXES.xlsx
# from https://population.un.org/wpp/Download/Standard/Population/
# These go 1.5 years later than Eurostat (mid 2020 cf start 2019) so projection not so important.

# Following https://www.actuaries.org.uk/system/files/field/document/CMI%20WP111%20v02%202019-04-25-%20Regular%20monitoring%20of%20England%20%20Wales%20population%20mortality.pdf
# and https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/articles/comparisonsofallcausemortalitybetweeneuropeancountriesandregions/januarytojune2020
# The former (actuaries) is more detailed, but using the notation of the latter (ONS), i.e., ASMR rather than SMR; csASMR etc.

# In this code, "ideal" day and week refer to idealised years that always have exactly 365 days and always start at the same week alignment (week 0 is 1-7 Jan).

import os,csv,sys,time,calendar,datetime
from collections import defaultdict
from subprocess import Popen,PIPE

# See https://ec.europa.eu/eurostat/cache/metadata/en/demomwk_esms.htm for country codes
# and https://population.un.org/wpp/Download/Metadata/Documentation/ for country names
countrycode='UK';countryname='United Kingdom'
#countrycode='FR';countryname='France'
#countrycode='ES';countryname='Spain'
#countrycode='CZ';countryname='Czechia'
#countrycode='IT';countryname='Italy'
#countrycode='SE';countryname='Sweden'
#countrycode='HU';countryname='Hungary'
meanyears=range(2015,2020)
agerange=(0,150)
targetyear=2020
# popsource = (Description, yearoffset, population [, projected population])
popsource=("WPP", 0.5, "WPP2019_POP_F15_1_ANNUAL_POPULATION_BY_AGE_BOTH_SEXES.csv")
#popsource=("Eurostat", 0.0, "urt_pjangrp3.tsv", "proj_19np.tsv")
useESP=False
update=True
mode="rASMR"# eASMR, rASMR, ecASMR or rcASMR

assert targetyear not in meanyears and 2020 not in meanyears
assert popsource[0] in ["WPP","Eurostat"]

deathsfn='demo_r_mwk_05.tsv'

print("Country:",countryname)
print("meanyears:",list(meanyears))
print("targetyear:",targetyear)
allyears=list(meanyears)+[targetyear]

# European Standard Population 2013
# https://webarchive.nationalarchives.gov.uk/20160106020035/http://www.ons.gov.uk/ons/guide-method/user-guidance/health-and-life-events/revised-european-standard-population-2013--2013-esp-/index.html
# ESP[i] = std pop for age group [5i,5(i+1)), except last one is [90,infinty)
ESP=[5000, 5500, 5500, 5500, 6000, 6000, 6500, 7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000, 4000, 2500, 1500, 1000]
nages=len(ESP)
agegrouprange=(agerange[0]//5,min((agerange[1]+4)//5,nages))# agegroup in [agegrouprage[0],agegrouprange[1])

# Convert Eurostat age string (Y_LT5, Y5-9, Y10-14, ..., Y85-89, Y_GE90) to integer 0-18
# Also convert Y_LT1, Y_1, Y_2, ..., Y_99, Y_GE100 to 0-18 according to which of the above brackets it belongs to
def age_s2i(s):
  assert s[0]=='Y'
  if s=='Y_LT1' or s=='Y_LT5': return 0
  i=0
  while i<len(s) and not s[i].isdigit(): i+=1
  j=i
  while j<len(s) and s[j].isdigit(): j+=1
  return min(int(s[i:j])//5,18)

# Inverse function to age_s2i
def age_i2s(i):
  if i==0: return 'Y_LT5'
  if i==nages-1: return 'Y_GE%d'%(5*i)
  return 'Y%d-%d'%(5*i,5*i+4)

# YYYY-MM-DD -> day number
def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

if update or not os.path.isfile(deathsfn):
  Popen("wget https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/demo_r_mwk_05.tsv.gz -O -|gunzip -c > %s"%deathsfn,shell=True).wait()

if popsource[0]=="Eurostat":
  for x in popsource[2:]:
    if not os.path.isfile(x):
      Popen("wget https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/"+x+".gz -O -|gunzip -c > "+x,shell=True).wait()
if popsource[0]=="WPP":
  if not os.path.isfile(popsource[2]):
    import pandas
    p=Popen("wget 'https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/"+noext(popsource[2])+".xlsx"+"' -O -",stdout=PIPE,shell=True)
    data=pandas.read_excel(p.stdout)
    data.to_csv(popsource[2], encoding='utf-8', index=False)
    p.wait()

# Convert ISO 8601 year,week to weighted list of year,day
# Resample leap years to range(365)
# (Not using datetime.date.fromisocalendar to avoid creating a dependancy on too recent a python version)
def isoweektodates(y,w):
  day=datetime.date.toordinal(datetime.date(y,1,1))+7*(w-2)
  while 1:
    dt=datetime.date.fromordinal(day)
    iso=datetime.date.isocalendar(dt)
    if iso[0]==y and iso[1]==w: break
    day+=1
  year2wd=defaultdict(int)# working days for each year present in week
  year2days=defaultdict(int)# total days for each year present in week
  for d in range(7):
    dt=datetime.date.fromordinal(day+d)
    year2days[dt.year]+=1
    # Only care about determining working days when crossing year boundary, in which case the week lies between 29 Dec and 6 Jan, inclusive.
    # Assume the only public holiday in this range is 1st Jan, or 2nd or 3rd Jan if 1st Jan falls on a Sunday, Saturday respectively.
    # Of course getting this wrong is only going to make a microscopic difference overall.
    if dt.weekday()>=6: holiday=True
    elif dt.month>1 or dt.day>3: holiday=False
    elif dt.day==1: holiday=True
    else: holiday=(dt.weekday()==0)
    if not holiday: year2wd[dt.year]+=1
  twd=sum(year2wd.values())
  l=defaultdict(float)
  for i in range(7):
    dt=datetime.date.fromordinal(day+i)
    wt=year2wd[dt.year]/twd/year2days[dt.year]
    day0=datetime.date.toordinal(datetime.date(dt.year,1,1))
    (y,d)=dt.year,day+i-day0
    # Resample leap year to range(365)
    if y%4:
      l[y,d]+=wt
    else:
      if d>0: l[y,d-1]+=d/366*wt
      if d<365: l[y,d]+=(365-d)/366*wt
  return [(y,d,l[y,d]) for (y,d) in l]

# Parse Eurostat death data
# End up with wd[age group number][year,idealday]=number of deaths for that age group, year, day
wd=[{} for age in range(nages)]
with open(deathsfn,'r') as fp:
  r=csv.reader(fp,delimiter='\t')
  first=True
  ok=[{y:set() for y in allyears} for age in range(nages)]
  for row in r:
    if first:
      assert row[0][:16]=='age,sex,unit,geo'
      # row[1:] is a list of strings like '2012W39 '; weeks are 01, 02, ..., 52, 53, 99.
      # Construct datelist and indexes
      dates=[]
      for i in range(len(row)-1,0,-1):
        y0,w0=map(int,row[i].strip().split('W'))
        if w0<=53:
          for (y,d,wt) in isoweektodates(y0,w0):
            if y in allyears:
              dates.append((i,y,d,wt))
      first=False
    else:
      (agestring,sex,unit,geo)=row[0].split(',')
      if sex=='T' and geo==countrycode and agestring[0]=='Y':
        age=age_s2i(agestring)
        for (i,y,d,wt) in dates:
          x=row[i].strip()
          if x[-1]=='p' or x[-1]=='e': x=x[:-1].strip()
          if x!=':':
            wd[age][y,d]=wd[age].get((y,d),0)+wt*float(x)
            ok[age][y].add(d)
  for age in range(nages):
    for y in meanyears:
      if len(ok[age][y])!=365: print("Missing data in year %d at age band %s"%(y,age_i2s(age)),file=sys.stderr);sys.exit(1)

lastidealday=min(len(ok[age][targetyear]) for age in range(nages))
numidealweeks=lastidealday//7# last+1 centered ideal week

def idealdaytostring(y,d0):
  d=int(d0*(365+(y%4==0))/365+.5)
  return datetime.date.fromordinal(datetime.date.toordinal(datetime.date(y,1,1))+d).strftime("%Y-%m-%d")

def int2(s):
  if s[-2:]==' p': return int(s[:-2])
  return int(s)

# End up with pp[year-minyear_pop][age group]=population for that year and age group, at proportion popsource[1] of the way through the year
if popsource[0]=="Eurostat":
  # Parse Eurostat population data
  with open(popsource[2],'r') as fp:
    r=csv.reader(fp,delimiter='\t')
    first=True
    for row in r:
      if first:
        assert row[0].strip()=='unit,sex,age,terrtypo,geo\\time'
        # Expecting row[1:] = 2019, 2018, 2017, 2016, 2015, 2014
        nyears_pop=len(row)-1
        minyear_pop=int(row[-1])
        pp=[[0]*nages for y in range(nyears_pop+2)]
        first=False
      else:
        (unit,sex,agestring,terr,geo)=row[0].split(',')
        # TOTAL territory not available, so construct it from URB+RUR+INT
        if sex=='T' and geo==countrycode and agestring[0]=='Y' and (terr=='URB' or terr=='RUR' or terr=='INT'):
          assert len(row)==nyears_pop+1
          assert ': z' not in row# Insist all data present
          age=age_s2i(agestring)
          for i in range(nyears_pop): pp[i][age]+=int2(row[nyears_pop-i])
  if countrycode=="UK":
    # Brexit means no Eurostat population projections for the UK
    pass
  else:
    # Parse Eurostat projected population data
    # Ages are Y_LT1, Y1, Y2, ..., Y99, Y_GE100, TOTAL
    with open(popsource[3],'r') as fp:
      r=csv.reader(fp,delimiter='\t')
      first=True
      for row in r:
        if first:
          assert row[0].strip()=='projection,unit,sex,age,geo\\time'
          # In 2020 expect row[-3:] = 2021, 2020, 2019
          for col in range(1,len(row)):
            if int(row[col].strip())==minyear_pop+nyears_pop: break
          else: raise LookupError("Year %d not found in %s"%(minyear_pop+nyears_pop,popsource[3]))
          first=False
        else:
          (proj,unit,sex,agestring,geo)=row[0].split(',')
          if proj=="BSL" and sex=='T' and geo==countrycode and agestring[0]=='Y':
            age=age_s2i(agestring)
            for i in range(2):
              pp[nyears_pop+i][age]+=int(row[col-i])
    nyears_pop+=2
else:
  # Parse UN WPP population data
  with open(popsource[2],'r') as fp:
    r=csv.reader(fp)
    start=True
    for row in r:
      if start:
        if row[0]=='Index':
          for (c,x) in enumerate(row):
            if x[:6]=='Region': cc=c
            if x[:9]=='Reference': yearcol=c
            if x=='0-4': agecol0=c
          agecol1=len(row)
          pp=[]
          start=False
      else:
        if row[cc]==countryname:
          y=int(row[yearcol])
          if len(pp)==0: minyear_pop=y
          assert y==minyear_pop+len(pp)# Insist years are contiguous
          l=[int(float(row[i].replace(' ',''))*1000+.5) for i in range(agecol0,agecol1)]
          l=l[:nages-1]+[sum(l[nages-1:])]# If the 90+ age group is subdivided then amalgamate it into one group
          pp.append(l)
  nyears_pop=len(pp)

# Number of deaths at age group a, year y0, in a 1 week interval centered around day d0 (0-364)
def D(y0,d0,a):
  t=0
  for x in range(y0*365+d0-3,y0*365+d0+4):
    y=x//365;d=x%365
    t+=wd[a][y,d]
  return t

# Estimated population at age group a, year y, ideal day d (0-364)
def E(y,d,a):
  yf=y+d/365-popsource[1]-minyear_pop# convert to ideal index of pp[]
  y0=min(int(yf),nyears_pop-2)
  y1=y0+1
  return (y1-yf)*pp[y0][a]+(yf-y0)*pp[y1][a]

# Print populations by age
if 1:
  print()
  print("                "+"        ".join("%4d"%y for y in allyears))
  s={y:sum(E(y,0,a) for a in range(nages)) for y in allyears}
  for a in range(nages): 
    print("%8s"%age_i2s(a),end="")
    #for y in allyears: n=E(y,0,a);print("   %9d"%n,end="")
    for y in allyears: n=E(y,0,a);print("   %8.2f%%"%(n/s[y]*100),end="")
    print()
  print("   TOTAL   ",end="")
  print("   ".join("%9d"%s[y] for y in allyears))
  print()

ASMR={y:[] for y in allyears}#      ASMR[y][w] = ASMR at year y, ideal week w (0-51)
cASMR={y:[0] for y in allyears}# cASMR[y][w+1] = cASMR at year y, ideal week w (0-51)
if useESP: REFPOP=ESP
for y in allyears:
  numw=numidealweeks if y==targetyear else 52
  for w in range(numw):
    POP=[E(2020,3+w*7,a) for a in range(nages)]
    if not useESP: REFPOP=POP
    d=3+w*7
    t=0
    for a in range(agegrouprange[0],agegrouprange[1]):
      t+=D(y,d,a)/E(y,d,a)*REFPOP[a]
    ASMR[y].append(t/sum(REFPOP)*sum(POP))
    cASMR[y].append(cASMR[y][-1]+ASMR[y][-1])

ASMR_bar=[]
for w in range(52):
  t=0
  for y in meanyears: t+=ASMR[y][w]
  ASMR_bar.append(t/len(meanyears))

cASMR_bar=[]
for w in range(53):
  t=0
  for y in meanyears: t+=cASMR[y][w]
  cASMR_bar.append(t/len(meanyears))

rcASMR=[]
ecASMR=[]
for w in range(numidealweeks+1):
  ecASMR.append(cASMR[targetyear][w]-cASMR_bar[w])
  rcASMR.append((cASMR[targetyear][w]-cASMR_bar[w])/cASMR_bar[52])

rASMR=[]
eASMR=[]
for w in range(numidealweeks):
  eASMR.append(ASMR[targetyear][w]-ASMR_bar[w])
  rASMR.append((ASMR[targetyear][w]-ASMR_bar[w])/(cASMR_bar[52]/52))

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

def noext(s):
  r=s.rfind('.')
  if r>=0: return s[:r]
  else: return s
  
def escape(s):
  return noext(s).replace('_','\\\_')

fn=countryname.replace(' ','')+'_'+mode+'.png'
po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
write('set terminal pngcairo font "sans,13" size 1920,1280')
write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
write('set output "%s"'%fn)
write('set key left')
title="Mortality in %s for %d"%(countryname,targetyear)
title+=' compared with %d-year average'%len(meanyears)+' for corresponding week of year\\n'
title+='Using '+("" if useESP else "dynamic variant of ")+mode+' measure.'
title+='   Age range %d - '%(agegrouprange[0]*5)+('%d.'%(agegrouprange[1]*5) if agegrouprange[1]<nages else 'infinity.')
title+='   Last date: %s\\n'%idealdaytostring(targetyear,lastidealday)
if popsource[0]=="Eurostat":
  title+='Sources: Eurostat '+', '.join(map(escape,popsource[2:]+(deathsfn,)))
else:
  title+='Sources: '+escape(popsource[2])+', Eurostat '+escape(deathsfn)
write('set title "%s"'%title)
write('set grid xtics lc rgb "#e0e0e0" lt 1')
write('set grid ytics lc rgb "#e0e0e0" lt 1')
write('set xtics nomirror')
write('set y2tics mirror')
write('set xtics rotate by 45 right offset 0.5,0')
write('set xdata time')
write('set format x "%Y-%m"')
write('set timefmt "%Y-%m-%d"')
if 'r' in mode:
  if 'c' in mode:
    graphlabel='excess deaths as % of expected deaths over year'
  else:
    graphlabel='excess deaths over week as % of expected deaths over week'
else:
  if 'c' in mode:
    graphlabel='cumulative excess deaths'
  else:
    graphlabel='excess deaths per week'
write('plot 0 lw 3 title "Baseline", "-" using 1:2 w lines lw 3 title "'+mode+' ('+graphlabel+')"')
if mode=="rASMR": source=[x*100 for x in rASMR]
if mode=="eASMR": source=eASMR
if mode=="rcASMR": source=[x*100 for x in rcASMR]
if mode=="ecASMR": source=ecASMR
for w in range(numidealweeks+('c' in mode)): write(idealdaytostring(targetyear,3*('c' in mode)+w*7),source[w])
write("e")
p.close()
po.wait()
print("Written %s"%fn)
