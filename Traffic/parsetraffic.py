from PIL import Image
import numpy as np
from os.path import join
from os import listdir

tdir='.'

# Google maps traffic levels
whitecol   = 0xffffff
greencol   = 0x63d668
orangecol  = 0xff974d
redcol     = 0xf23c32
darkredcol = 0x811f1f

levels={}
cachefn=join(tdir,'trafficlevels')
try:
  with open(cachefn,'r') as fp:
    for x in fp:
      if x[:1]!='#':
        y=x.strip().split()
        levels[(y[0],y[1])]=[int(z) for z in y[2:]]
except FileNotFoundError:
  pass

for fn in sorted(listdir(tdir)):
  if fn[:6]=='London' and fn[-4:]=='.png':
    f=fn[:-4].find('.')
    key=(fn[:f],fn[f+1:-4])
    if key not in levels:
      im_frame = Image.open(join(tdir,fn))
      np_frame = np.array(im_frame)
      im_frame.close()
      rgb = np_frame[:,:,0]*65536+np_frame[:,:,1]*256+np_frame[:,:,2]
      numcols = [np.sum(rgb==v) for v in [whitecol,greencol,orangecol,redcol,darkredcol]]
      if numcols[0]<1000000:
        levels[key]=numcols
        print("Read",fn)
      else:
        print("Discarding",fn)

# Work out overall congestion as weighted sum of orange, red and dark red
ll=sorted(levels);n=len(ll)
l0=[]
for x in sorted(ll):
  m=levels[x]
  c=m[2]+3*m[3]+10*m[4]
  l0.append(c)

# Work out cumulative version of l0 with padded ends
l0pad=l0[:168]+l0+l0[-168:]
cum=[]
t=0
for c in l0pad:
  cum.append(t)
  t+=c
cum.append(t)

# Create 24-hour and 168-hour moving averages
l1=[];l2=[]
for (a,l) in [(24,l1),(168,l2)]:
  assert a%2==0;a2=a//2
  for i in range(168,n+168):
    l.append(int((cum[i+a2]-cum[i-a2+1]+(l0pad[i-a2]+l0pad[i+a2])/2+a2)/a))

with open(cachefn,'w') as fp:
  print("# Location      LocalDateTime    White    Green   Orange      Red  Darkred     Cong   Cong24  Cong168",file=fp)
  for i in range(n):
    (loc,date)=ll[i]
    print(loc,date,' '.join("%8d"%x for x in levels[(loc,date)][:5]),"%8d %8d %8d"%(l0[i],l1[i],l2[i]),file=fp)
