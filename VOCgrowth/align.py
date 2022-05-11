from __future__ import division, print_function
# Aligns SARS-CoV-2 fasta files to reference genome Wuhan-Hu-1
# Use pypy for speed

import sys
from math import sqrt

with open('refgenome') as fp:
  ref=fp.read().strip()

N=len(ref)
r=9
refdict={}
for i in range(N-r+1): refdict.setdefault(ref[i:i+r],[]).append(i)

fp=sys.stdin
try:
  name=next(fp).rstrip()
  last=False
except StopIteration:
  last=True

while not last:
  lines=[]
  while 1:
    try:
      line=next(fp).rstrip()
      if line[:1]=='>': name=line;break
      lines.append(line)
    except StopIteration:
      last=True
      break
  genome=''.join(lines)
  
  # Offset = (index in ref genome) - (index in current genome)
  M=len(genome)
  prev=0
  offsetcount={}
  for i in range(M-r+1):
    #print(i,refdict.get(genome[i:i+r],[]))
    k=genome[i:i+r]
    if k in refdict:
      for x in refdict[k]:
        offsetcount[x-i]=offsetcount.get(x-i,0)+1
  
  okoffsets=set(x for x in offsetcount if offsetcount[x]>=20)
  noff=len(okoffsets)
  
  #for x in okoffsets:
  #  print(x,offsetcount[x])
  
  infinity=65535
  offsets=[[infinity,infinity] for i in range(M)]
  
  # Approach from right
  nearest=infinity
  for i in range(M-r,-1,-1):
    k=genome[i:i+r]
    if k in refdict:
      for x in refdict[k]:
        if x-i in okoffsets: nearest=x-i
    offsets[i][1]=nearest
  
  # Approach from left
  nearest=infinity
  for i in range(r-1,M):
    k=genome[i-r+1:i+1]
    if k in refdict:
      for x in refdict[k]:
        if x-(i-r+1) in okoffsets: nearest=x-(i-r+1)
    offsets[i][0]=nearest
  
  st=[[infinity,0],[infinity,0]]
  bp=[[0,0] for i in range(M)]
  for i in range(M):
    newst=[[infinity,1000000000],[infinity,1000000000]]
    for j in range(2):
      off=offsets[i][j]
      best=[1000000000,]
      if off!=infinity:
        for k in range(2):
          [p,v]=st[k]
          v1=v+int(sqrt(abs(off-p)))*(p!=infinity)+(genome[i]!=ref[i+off])
          if v1<best[0]: best=[v1,k]
        bp[i][j]=best[1]
        newst[j]=[off,best[0]]
    st=newst
    if 0:
      print(i,genome[i],end='   ')
      for j in range(2):
        off=offsets[i][j]
        print(off,end=' ')
        if off!=infinity: print(ref[i+off]+'*'*(ref[i+off]==genome[i]),bp[i][j],end='   ')
        else: print('-',end='   ')
      print(st)
  
  best=[1000000000,]
  for k in range(2):
    [p,v]=st[k]
    if v<best[0]: best=[v,k]
  s=best[1]
  offs=[infinity]*M
  for i in range(M-1,-1,-1):
    offs[i]=offsets[i][s]
    s=bp[i][s]
    
  if 0:
    for i in range(M):
      off=offs[i]
      print("%5d"%i,genome[i],"   ",off,end=' ')
      if off!=infinity: print(ref[i+off]+'!'*(genome[i] in 'ACGT' and ref[i+off]!=genome[i]))
      else: print('-')
  
  print(name)
  out=['N']*N
  for i in range(M):
    off=offs[i]
    if off!=infinity: out[i+off]=genome[i]
  print(''.join(out))
