# Simple Zoe England-wide symptom-based prevalence from map data

import json,os,sys
import numpy as np
from stuff import *

mindate='2020-10-01'
#mindate='2021-03-01'
country='England'

print('Country:',country)

np.set_printoptions(precision=3,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=100000)

# Convert Zoe into convenient sequence
# zvals[day]
minday=datetoday(mindate)
maxt=-1
zvals=np.zeros(10000)
zdates=np.zeros(10000,dtype=bool)
for date in os.listdir('zoemapdata'):
  if date[:2]!='20': continue# Ignore spurious files
  t=datetoday(date)-minday
  if t>=0:
    maxt=max(maxt,t)
    with open(os.path.join('zoemapdata',date),'r') as fp:
      zm=json.load(fp)
      zdates[t]=1
      for d in zm.values():# Loop over locations
        if (d['country'] if 'country' in d else d['country_x'])==country: zvals[t]+=d['predicted_covid_positive_count']/d['respondent']*d['population']
n=maxt+1
# Interpolate missing dates
for t in range(n):
  if not zdates[t]:
    for t0 in range(t-1,-1,-1):
      if zdates[t0]: break
    else: raise RuntimeError("Missing Zoe entry at start")
    for t1 in range(t+1,n):
      if zdates[t1]: break
    else: raise RuntimeError("Missing Zoe entry at end")
    zvals[t]=((t1-t)*zvals[t0]+(t-t0)*zvals[t1])/(t1-t0)

with open('zoeprev','w') as fp:
  for t in range(n):
    print(daytodate(minday+t),zvals[t],file=fp)

print('gnuplot')
print('set timefmt "%Y-%m-%d";set format x "%Y-%m-%d";set xdata time')
print('plot "zoeprev" u 1:2 w linespoints lw 3 title "Zoe symptom-based prevalence in England"')
