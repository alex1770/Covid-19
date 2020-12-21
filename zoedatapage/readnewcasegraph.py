# Conversion of Zoe new case graphs to numbers

# Input (eg):
# London.2020-11-19.ppm 98  1027 73  420     196  958   96  420  140    0   2020-06-01  2020-11-01
# London.2020-12-19.ppm 99  1028 74  418     182  967   95  420  140    0   2020-06-01  2020-12-01
# filename              x0  x1   y0   y1     x2   x3    y2   y3   v2   v3   date_x2     date_x3
#           top-leftian         top  bot               top  bot
#                       ------------------   ------------------  ---------------------------------
#                       Region of interest   Calibration points  Calibration values

# Disc:
# Colour 1f77b4
# Diameter 12ish (11 inner + 1 partial on each side)
# Favour discs with lots of this col within rad0 and as little as poss between rad1 and rad2
col=[0x1f,0x77,0xb4]
rad0,rad1,rad2=6,7,8

import sys,time,calendar
import numpy as np
from PIL import Image

def datetoday(x):
  t=time.strptime(x+'UTC','%Y-%m-%d%Z')
  return calendar.timegm(t)//86400

def daytodate(r):
  t=time.gmtime(r*86400)
  return time.strftime('%Y-%m-%d',t)

[x0,x1,y0,y1,x2,x3,y2,y3,v2,v3]=[float(x) for x in sys.argv[2:12]]
date_x2=sys.argv[12]
date_x3=sys.argv[13]

imf=Image.open(sys.argv[1])
npf=np.array(imf)
imf.close()
hei,wid,ncol=npf.shape
isin=np.linalg.norm((npf.astype(int)-col),axis=2)<10

d2=datetoday(date_x2)
d3=datetoday(date_x3)

def pixtoday(x): return (x-x2)/(x3-x2)*(d3-d2)+d2
def daytopix(d): return (d-d2)/(d3-d2)*(x3-x2)+x2
def pixtoval(y): return (y-y2)/(y3-y2)*(v3-v2)+v2

d0=int(pixtoday(x0)+.999)
d1=int(pixtoday(x1))

for day in range(d0,d1+1):
  date=daytodate(day+4)# To make it into publication date
  x=int(daytopix(day)+.5)
  best=(-1e9,)
  for y in range(int(y0),int(y1)+1):
    sc=0
    for dx in range(-rad2,rad2+1):
      for dy in range(-rad2,rad2+1):
        r2=dx**2+dy**2
        if r2<rad0**2: sc+=isin[y+dy,x+dx]*(1+(dy<0))
        elif r2>=rad1**2 and r2<rad2**2: sc-=isin[y+dy,x+dx]
    #print(y,sc)
    if sc>best[0]: best=(sc,y)
  y=best[1]
  #print(date,x,y,"%.1f"%(pixtoval(y)*10))
  print(date,"%6.1f"%(pixtoval(y)*10))# Turn it into cases per million
  npf[y,x]=[255,255,255]
  if day%10==0:
    imf=Image.fromarray(npf)
    imf.save('temp.png')

imf=Image.fromarray(npf)
imf.save('temp.png')
