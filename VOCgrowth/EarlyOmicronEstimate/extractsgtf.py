from stuff import *

# Get weekday pattern from case data in order to identify exact date on SGTF graph
# 0 mod 7 is Thursday in daytodate notation (being 1970-01-01)
nc={}
with open('SAcases','r') as fp:
  for x in fp:
    y=x.split()
    nc[datetoday(y[0])]=int(y[1])

minday=min(nc)
maxday=max(nc)
c0=[0]*7
c1=[0]*7
for d in range(minday+3,maxday-3):
  ex=[nc[r] for r in range(d-3,d+4)]
  if min(ex)>=50:
    i=d%7
    c0[i]+=1
    c1[i]+=nc[d]*7/sum(ex)

#for i in range(7):
#  print(i,c1[i]/c0[i])
# Thur 1.184
# Fri  1.170
# Sat  1.122
# Sun  0.913
# Mon  0.655
# Tue  0.766
# Wed  1.158

if 0:
  infile='OmicronSGTF.png'
  dateorigin=datetoday('2021-10-01')-564
  row0,row1=23,359
  col0,col1=81,614
  y0=(0,358);y1=(50,43)
  z0=(0,357);z1=(1600,126)

if 1:
  infile='OmicronSGTF_frompdf.png'
  dateorigin=datetoday('2021-10-01')-564
  row0,row1=11,345
  col0,col1=81,614
  y0=(0,344.5);y1=(50,32)
  z0=(0,344.5);z1=(2000,57.5)

# SGTF image from slide 12 of https://sacoronavirus.co.za/2021/11/25/sars-cov-2-sequencing-new-variant-update-25-november-2021/
# resized down by a factor of 2/3 in order to get 1 horizontal pixel = 1 day.
from PIL import Image
import numpy as np
im_frame = Image.open(infile)
cc = np.array(im_frame,dtype=int)
im_frame.close()
# Top-leftian, row before column

r=cc.shape[0]
c=cc.shape[1]

# Get blueness
bb=cc[:,:,2]*2-(cc[:,:,0]+cc[:,:,1])

def process(bb,name):
  bb1=bb[row0:row1,:]
  mm=row0+np.argmax(bb1,axis=0)
  im=Image.fromarray(((bb-bb.min())/(bb.max()-bb.min())*255.999+0.0005).astype(np.dtype('uint8')))
  im.save(name+'_filtered.png')

  oo=cc.astype(np.dtype('uint8'))
  for x in range(col0,col1): oo[mm[x],x]=[255,0,0]
  im=Image.fromarray(oo)
  im.save(name+'_sgtf.png')

  sgtf={}
  for x in range(col0,col1):
    sgtf[daytodate(dateorigin+x)]=(mm[x]-y1[1])/(y0[1]-y1[1])*(y0[0]-y1[0])+y1[0]
  with open(name+'_sgtf','w') as fp:
    for date in sorted(list(sgtf)):
      print(date,"%6.2f"%sgtf[date],file=fp)

  return mm,sgtf

process(bb,'simple')

lrantialias=bb-np.maximum(np.roll(bb,1,1),np.roll(bb,-1,1))
process(lrantialias,'LRantialias')

# Hybrid because deantialiasing method is likely to work well for the vertical spike, but not when derivative is low.
spike=605
hybrid=np.concatenate([bb[:,:spike],lrantialias[:,spike:]],axis=1)
mm,sgtf=process(hybrid,'hybrid')

dd=cc[:,:,0]-np.maximum(cc[:,:,1],cc[:,:,2])
oo=(dd>3).astype(np.dtype('uint8'))*255
im=Image.fromarray(oo)
im.save('temp.png')

ee=(dd>3)*1000+np.tile(np.arange(r-1,-1,-1)[:,None],(1,c))
process(ee,'simplered')

oo=cc.astype(np.dtype('uint8'))
nn=np.zeros(c)
for x in range(col0,col1):
  s0=1
  s1=10
  f=0.5
  mx=0
  for y in range(row1-1,row0-1,-1):
    if abs(y-mm[x])>1:
      s0=(1-f)*s0+f*1
      s1=(1-f)*s1+f*dd[y,x]
      #print(y,dd[y,x],s1/s0)
      if s1/s0>5: mx=y
  nn[x]=mx
  oo[mx,x]=[0,255,0]
  oo[mm[x],x]=[255,0,0]

im=Image.fromarray(oo)
im.save('sgtf+counts.png')

with open('SA_sgtf','w') as fp:
  print("#     Date  %SGTF   Tests  num(S-)  num(S+)",file=fp)
  for x in range(col0,col1):
    if nn[x]>0:
      date=daytodate(dateorigin+x)
      n=max((nn[x]-z1[1])/(z0[1]-z1[1])*(z0[0]-z1[0])+z1[0],0)
      s=sgtf[date]
      print(date,"%6.2f  %6.1f   %6.1f   %6.1f"%(s,n,s/100*n,(1-s/100)*n),file=fp)
