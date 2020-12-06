# From https://www.gov.uk/guidance/the-r-number-in-the-uk#history
a="""
The R number range for the UK is 0.8-1.0 and the growth rate range is -3% to -1% per day as of 4 December 2020.
The R number range for the UK is 0.9-1.0 and the growth rate range is -2% to 0% per day as of 27 November 2020.
The R number range for the UK is 1.0-1.1 and the growth rate range is 0% to +2% per day as of 20 November 2020.
The R number range for the UK is 1.0-1.2 and the growth rate range is +1% to +3% per day as of 13 November 2020.
The R number range for the UK is 1.1-1.3 and the growth rate range is +2% to +4% per day as of 6 November 2020.
The R number range for the UK is 1.1-1.3 and the growth rate range is +2% to +4% per day as of 30 October 2020.
The R number range for the UK is 1.2-1.4 and the growth rate range is +3% to +6% per day as of 23 October 2020.
The R number range for the UK is 1.3-1.5 and the growth rate range is +4% to +7% per day as of 16 October 2020.
The R number range for the UK is 1.2-1.5 and the growth rate range is +4% to +9% per day as of 9 October 2020.
The R number range for the UK is 1.3-1.6 and the growth rate range is +5% to +9% per day as of 2 October 2020.
The R number range for the UK is 1.2-1.5 and the growth rate range is +4% to +8% per day as of 25 September 2020.
The R number range for the UK is 1.1-1.4 and the growth rate range is +2% to +7% per day as of 18 September 2020.
The R number range for the UK is 1.0-1.2 and the growth rate range is -1% to +3% per day as of 11 September 2020.
The R number range for the UK is 0.9-1.1 and the growth rate range is -1% to +2% per day as of 4 September 2020.
The R number range for the UK is 0.9-1.1 and the growth rate range is -2% to +1% per day as of 28 August 2020.
The R number range for the UK is 0.9-1.1 and the growth rate range is -3% to +1% per day as of 21 August 2020.
The R number range for the UK is 0.8-1.0 and the growth rate range is -4% to -1% per day as of 14 August 2020.
The R number range for the UK is 0.8-1.0 and the growth rate range is 0% to -5% per day as of 7 August 2020.
The R number range for the UK is 0.8-0.9 and the growth rate range is -1% to -4% as of 31 July 2020.
The R number range for the UK is 0.7-0.9 and the growth rate range is -4% to -1% as of 24 July 2020.
The R number range for the UK is 0.7-0.9 and the growth rate range is -5% to -1% as of 17 July 2020.
The R number range for the UK is 0.7-0.9 and the growth rate range is -5% to -2% as of 10 July 2020.
The R number range for the UK is 0.7-0.9 and the growth rate range is -6% to -0% as of 3 July 2020.
The R number range for the UK is 0.7-0.9 and the growth rate range is -4% to -2% as of 25 June 2020.
The R number range for the UK is 0.7-0.9 and the growth rate range is -4% to -2% as of 19 June 2020.
The R number range for the UK is 0.7-0.9 as of 12 June 2020.
The R number range for the UK is 0.7-0.9 as of 5 June 2020.
The R number range for the UK is 0.7-0.9 as of 29 May 2020.
The R number range for the UK is 0.7-1.0 as of 22 May 2020.
""".strip().split('\n')

from math import log,sqrt
import numpy as np
from scipy.optimize import minimize

# Fit alpha, beta to R = (1+lambda/beta)^alpha, where lambda=growh rate
LL=[];RR=[]
for x in a:
  y=x.split()
  if 'growth' in y:
    i=y.index('growth')
    ll=[]
    for j in [i+4,i+6]:
      ll.append(log(1+float(y[j][:-1])/100))
    ll.sort()# Because they've messed up and sometimes present higher before lower
    for l in ll: LL.append(l)
    for z in y[8].split('-'):
      RR.append(log(float(z)))

def do(beta):
  sxy=sxx=0
  for (l,r) in zip(LL,RR):
    x=log(1+l/beta)
    y=r
    sxx+=x*x;sxy+=x*y
  al=sxy/sxx;err=0
  for (l,r) in zip(LL,RR):
    err+=(r-al*log(1+l/beta))**2
  return err,al

res=minimize(lambda x:do(x)[0],[1],method='BFGS')
assert res.success
be=res.x[0]
al=do(be)[1]
#print(al,be)

offset=-30
import time
prev=None
for x in a[::-1]:
  y=x.split()
  rr=[float(z) for z in y[8].split('-')]
  dt=time.strftime('%Y-%m-%d',time.gmtime(time.mktime(time.strptime(y[-3]+' '+y[-2]+' '+y[-1][:-1]+' 12','%d %B %Y %H'))+offset*86400))
  R=sum(rr)/2
  if prev: print(dt+", %5.2f"%prev)
  print(dt+", %5.2f"%R)
  prev=R
