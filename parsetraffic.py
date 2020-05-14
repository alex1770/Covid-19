from PIL import Image
import numpy as np
from os.path import join
from os import listdir

tdir='Traffic'

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
        levels[y[0]]=[int(z) for z in y[1:]]
except FileNotFoundError:
  pass

for fn in sorted(listdir(tdir)):
  if fn[:6]=='London' and fn[-4:]=='.png' and fn[:-4] not in levels:
    im_frame = Image.open(join(tdir,fn))
    np_frame = np.array(im_frame)
    im_frame.close()
    rgb = np_frame[:,:,0]*65536+np_frame[:,:,1]*256+np_frame[:,:,2]
    numcols = [np.sum(rgb==v) for v in [whitecol,greencol,orangecol,redcol,darkredcol]]
    levels[fn[:-4]]=numcols
    print("Read",fn)
    #traffic = numcols[2]+3*numcols[3]+10*numcols[4]
    #print(fn,"%8d"%traffic)

with open(cachefn,'w') as fp:
  for name in sorted(levels):
    print(name,' '.join("%8d"%x for x in levels[name]),file=fp)

