# Quick and dirty conversion of Zoe incidence files to regional new case numbers (ranges)

import os
from math import sqrt

def mid(a,b):
  #return sqrt(a*b)
  return a**.55*b**.45

for fn in sorted(os.listdir('.')):
  if fn[:10]!='incidence.': continue
  with open(fn,'r') as fp:
    for x in fp:
      if x[:1].isdigit() and '-' in x:
        y=x.split('-')
        if len(y)==2: ran=(int(y[0]),int(y[1]))
      m='<span class="region_name">'
      n=len(m)
      if x[:n]==m:
        f=x.find('<',n);assert f>n
        reg=x[n:f]
        print(fn[10:20],"%7.1f "%mid(*ran),reg,"%5d %5d"%ran)
