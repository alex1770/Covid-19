from stuff import *
import sys
import numpy as np

def myfloat(x):
  try:
    return float(x)
  except:
    return 0

cols=[int(x) for x in sys.argv[1:]]
n=len(cols)

data={}
minday=1000000
maxday=-1
for x in sys.stdin:
  y=x.split()
  day=datetoday(y[0])
  data[day]=np.array([myfloat(y[c-1]) for c in cols])
  if day<minday: minday=day
  if day>maxday: maxday=day

z=np.zeros(n)
for day in range(minday,maxday+1):
  ave=sum(data.get(d,z) for d in range(day-3,day+4))/7
  print(daytodate(day),*['  %10.2f'%ave[i] for i in range(n)])
  
  
