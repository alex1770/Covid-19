#

from math import exp,log
import getdata,time

data=getdata.getcountrydata("UK")
dates=data[0]
cases=data[5]
deaths=data[6]
for i in range(len(cases)-1,-1,-1):
  if cases[i]<10: break
dates=dates[i+1:]
cases=cases[i+1:]
deaths=deaths[i+1:]
n=len(dates)

m=[]
for i in range(n-6):
  v=sum(cases[i+j] for j in range(7))/7
  m.append((dates[i],v))

n=len(m)
d=14
offset=3
#mgt=5.5
#alpha=2.29;beta=alpha/6.29
alpha=1.991;beta=0.336
for i in range(n-d):
  #print(m[i][0],m[i][1])
  dt=time.strftime('%Y-%m-%d',time.gmtime(time.mktime(time.strptime(m[i][0]+' 12','%Y-%m-%d %H'))+offset*86400))
  lam=log(m[i+d][1]/m[i][1])/d# growth rate
  #R=exp(lam*mgt)
  R=(1+lam/beta)**alpha
  print(dt+',',R)
