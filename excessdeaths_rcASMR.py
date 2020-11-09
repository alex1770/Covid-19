# Trying to replicate rcASMR calculation from bottom of https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/articles/comparisonsofallcausemortalitybetweeneuropeancountriesandregions/januarytojune2020
# Simplified version not handling week 53, leap years, or population projections to 2019, 2020 properly

# Population source: https://ec.europa.eu/eurostat/web/products-datasets/product?code=urt_pjangrp3 giving urt_pjangrp3.tsv
# (https://ec.europa.eu/eurostat/en/web/products-datasets/-/demo_pjangroup goes back further but doesn't have a 85-89 age group)
# For population projections (not currently used), can use
# EU: https://ec.europa.eu/eurostat/web/products-datasets/product?code=proj_19np
# UK: (uses midpoints of years): https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/tablea21principalprojectionukpopulationinagegroups
#     or                         https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/tablel21oldagestructurevariantukpopulationinagegroups
# but for simplicity not using these to start with (instead linearly interpolate off the end of the actuals)
# Going to be simpler to switch to using:
# https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2019_POP_F15_1_ANNUAL_POPULATION_BY_AGE_BOTH_SEXES.xlsx
# from https://population.un.org/wpp/Download/Standard/Population/
# (though I think the linear interpolation off the end isn't actually that bad an estimate)

# Mortality source: https://data.europa.eu/euodp/en/data/dataset/QrtzdXsI5w26vnr54SIzpQ giving demo_r_mwk_05.tsv
# Explanation: https://ec.europa.eu/eurostat/cache/metadata/en/demomwk_esms.htm
# Ignoring week 99, which appears to contain no data in the examples I've looked at.
# The week number (1-53) is the ISO 8601 week: weeks run Mon-Sun and are assigned to the year in which Thursday falls.
# For the moment just ignore week 53 and (slightly wrongly) pretend that week i is aligned at the same point in the year for any year.

import os,csv,time,calendar
from subprocess import Popen,PIPE

# See https://ec.europa.eu/eurostat/cache/metadata/en/demomwk_esms.htm for country codes
countrycode='UK';countryname='UK'
#countrycode='FR';countryname='France'
#countrycode='ES';countryname='Spain'
meanyears=range(2015,2020)
targetyear=2020
assert targetyear not in meanyears and 2020 not in meanyears
update=False

print("Country:",countryname)
print("meanyears:",list(meanyears))
print("targetyear:",targetyear)
allyears=list(meanyears)+[targetyear]
minyear=min(allyears)

# YYYY-MM-DD -> day number
def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

deathsfn='demo_r_mwk_05.tsv'
if update or not os.path.isfile(deathsfn):
  if os.path.exists(deathsfn): os.remove(deathsfn)
  Popen("wget https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/demo_r_mwk_05.tsv.gz -O -|gunzip -c > %s"%deathsfn,shell=True).wait()

popfn='urt_pjangrp3.tsv'
if update or not os.path.isfile(popfn):
  if os.path.exists(popfn): os.remove(popfn)
  Popen("wget https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/urt_pjangrp3.tsv.gz -O -|gunzip -c > %s"%popfn,shell=True).wait()

# Ages: Y_LT5, Y5-9, Y10-14, Y15-19, Y20-24, Y25-29, Y30-34, Y35-39, Y40-44, Y45-49, Y50-54, Y55-59, Y60-64, Y65-69, Y70-74, Y75-79, Y80-84, Y85-89, Y_GE90
# (Also TOTAL, UNK)

wd={}
minr=1000000;maxr=0
with open(deathsfn,'r') as fp:
  r=csv.reader(fp,delimiter='\t')
  first=True
  for row in r:
    if first:
      assert row[0][:16]=='age,sex,unit,geo'
      # row[1:] is a list of strings like '2012W39 '; weeks are 01, 02, ..., 52, 53, 99.
      # Construct datelist and indexes
      dates=[]
      for i in range(len(row)-1,0,-1):
        y,w=map(int,row[i].strip().split('W'))
        if y>=minyear and w<53:# Just ignore week 53 for the moment. Perhaps treat it properly later.
          dates.append((i,y,w,daytodate(datetoday("%4d-01-07"%y)+7*(w-1))))
      first=False
    else:
      (age,sex,unit,geo)=row[0].split(',')
      if sex=='T' and geo==countrycode and age[0]=='Y':
        if age not in wd: wd[age]=[-1]*len(dates)
        for (r,(i,y,w,d)) in enumerate(dates):
          x=row[i].strip()
          if x[-1]=='p' or x[-1]=='e': x=x[:-1].strip()
          if x!=':':
            wd[age][r]=int(x)
            if r<minr: minr=r
            if r>=maxr: maxr=r+1


# Reduce to valid date range and check contiguous
dd={}
for age in wd:
  dd[age]=wd[age][minr:maxr]
  assert -1 not in dd[age]
dates=dates[minr:maxr]
N=maxr-minr

def int2(s):
  if s[-2:]==' p': return int(s[:-2])
  return int(s)

pp={}
with open(popfn,'r') as fp:
  r=csv.reader(fp,delimiter='\t')
  first=True
  for row in r:
    if first:
      assert row[0]=='unit,sex,age,terrtypo,geo\\time'
      # Expecting row[1:] = 2019, 2018, 2017, 2016, 2015, 2014
      nyears=len(row)-1
      minyear=int(row[-1])
      # Something
      first=False
    else:
      (unit,sex,age,terr,geo)=row[0].split(',')
      # TOTAL territory not available, so construct it from URB+RUR+INT
      if sex=='T' and geo==countrycode and age[0]=='Y'and (terr=='URB' or terr=='RUR' or terr=='INT'):
        if age not in pp: pp[age]=[0]*nyears
        assert len(row)==nyears+1
        assert ': z' not in row# Insist all data present
        for i in range(nyears): pp[age][i]+=int2(row[nyears-i])

# Number of deaths at year y, week k (1-52), age group a
def D(y,w,a):
  return dd[a][(y-dates[0][1])*52+w-1]
    
# Estimated population at year y, week k, age group a
def E(y,w,a):
  #return ESP[age_s2i(a)]/sum(ESP)*67e6#test
  yy=y-minyear
  if yy<nyears-1: y0=yy;y1=yy+1
  else: y0=nyears-2;y1=nyears-1
  yf=yy+(w-.5)/52
  return (y1-yf)*pp[a][y0]+(yf-y0)*pp[a][y1]

# European Standard Population 2013
# https://webarchive.nationalarchives.gov.uk/20160106020035/http://www.ons.gov.uk/ons/guide-method/user-guidance/health-and-life-events/revised-european-standard-population-2013--2013-esp-/index.html
# ESP[i] = std pop for age group [5i,5(i+1)), except last one is [90,infinty)
ESP=[5000, 5500, 5500, 5500, 6000, 6000, 6500, 7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000, 4000, 2500, 1500, 1000]
nages=len(ESP)
    
def age_s2i(s):
  assert s[0]=='Y'
  if s=='Y_LT5': return 0
  i=0
  while i<len(s) and not s[i].isdigit(): i+=1
  j=i
  while j<len(s) and s[j].isdigit(): j+=1
  return int(s[i:j])//5

def age_i2s(i):
  if i==0: return 'Y_LT5'
  if i==nages-1: return 'Y_GE%d'%(5*i)
  return 'Y%d-%d'%(5*i,5*i+4)

# Print populations by age
if 1:
  print()
  print("                "+"        ".join("%4d"%y for y in allyears))
  s={y:0 for y in allyears}
  for a in range(nages): 
    aa=age_i2s(a)
    print("%8s"%aa,end="")
    for y in allyears: n=E(y,0.5,aa);s[y]+=n;print("   %9d"%n,end="")
    print()
  print("   TOTAL   ",end="")
  print("   ".join("%9d"%s[y] for y in allyears))
  print()

ASMR={}#  ASMR[y][w-1] = ASMR at year y, week w (1-52)
cASMR={}# cASMR[y][w] = cASMR at year y, week w (1-52)
for y in allyears:
  if y<2020: numw=52
  else: numw=dates[N-1][2]
  ASMR[y]=[]
  cASMR[y]=[0]
  #print(y,end="")
  td=[0]*nages;te=[0]*nages
  for w in range(1,numw+1):
    t=0
    for a in range(nages):
      aa=age_i2s(a)
      t+=D(y,w,aa)/E(y,w,aa)*ESP[a]#E(2020,40,aa)
      #if w==8: print(" %8.0f"%(t/sum(ESP)*1e6),end="")
      #if w==8: print(" %8.4f"%(D(y,w,aa)/E(y,w,aa)*100),end="")
      td[a]+=D(y,w,aa)
      te[a]+=E(y,w,aa)
    #if w==8: print()
    #if w==17: print(''.join(" %8.4f"%(100*d/e) for (d,e) in zip(td,te)))
    ASMR[y].append(t/sum(ESP))
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
numw=dates[N-1][2]
for w in range(numw+1):
  rcASMR.append((cASMR[targetyear][w]-cASMR_bar[w])/cASMR_bar[52])

dASMR=[]
for w in range(numw):
  dASMR.append(ASMR[targetyear][w]-ASMR_bar[w])

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

fn=countryname+'_rcASMR.png'
po=Popen("gnuplot",shell=True,stdin=PIPE);p=po.stdin
write('set terminal pngcairo font "sans,13" size 1920,1280')
write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
write('set output "%s"'%fn)
write('set key left')
#write('set logscale y')
title="Mortality in %s for %d"%(countryname,targetyear)
title+=' compared with %d-year average'%len(meanyears)+' for corresponding week of year, using rcASMR measure\\n'
title+='Last date: %s. '%(dates[-1][3])
title+='Sources: Eurostat urt\\\_pjangrp3 and demo\\\_r\\\_mwk\\\_05'
write('set title "%s"'%title)
write('set grid xtics lc rgb "#e0e0e0" lt 1')
write('set grid ytics lc rgb "#e0e0e0" lt 1')
write('set xtics nomirror')
write('set y2tics mirror')
write('set xtics rotate by 45 right offset 0.5,0')
write('set xdata time')
write('set format x "%Y-%m"')
write('set timefmt "%Y-%m-%d"')
write('plot "-" using 1:2 w lines title "rcASMR"')
#for w in range(numw): write(dates[(2020-dates[0][1])*52+w][3],dASMR[w]*sum(ESP))
for w in range(numw+1): write(dates[(2020-dates[0][1])*52+w-1][3],rcASMR[w]*100)
write("e")
p.close()
po.wait()
print("Written %s"%fn)

#for y in allyears:
#  print(y,"%.10f  %.10f"%(ASMR[y][8],ASMR[y][8]/ASMR_bar[8]))
