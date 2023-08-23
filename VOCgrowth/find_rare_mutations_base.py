# Input = aligned fasta, no linebreaks in a genome
# Program will look for positions with good smoothed likelihood ratios:
# E.g., ((number in example set with C in position 21612) / (number in example set without C in position 21612)) / 
#       ((number in all genomes with C in position 21612) / (number in all genomes without C in position 21612))
# where "with C" means "could have C", i.e., not {A, G or T}
# and "without C" means "couldn't have C", i.e., A or G or T.
# I.e., non-{A,C,G,T} act as wildcards

import sys
from stuff import getaa
from math import log

allcounts=[]
with open("nucleotide_counts") as fp:
  for x in fp:
    y=x.split()
    d={}
    for z in y[1:]:
      d[z[0]]=int(z[1:])
    allcounts.append(d)
n=len(allcounts)
alltot=[sum(x.values()) for x in allcounts]
allwild=[sum(x[y] for y in x if y not in "ACGT") for x in allcounts]

excounts=[{} for i in range(n)]
header=[]
example=[]
for x in sys.stdin:
  x=x.rstrip()
  if x[0]=='>': header.append(x);continue
  example.append(x)
  for (i,c) in enumerate(x):
    excounts[i][c]=excounts[i].get(c,0)+1
extot=[sum(x.values()) for x in excounts]
exwild=[sum(x[y] for y in x if y not in "ACGT") for x in excounts]

#disc=len([x for x in intersection if x=="#"])
#print(disc,"discrepanc"+["ies","y"][disc==1])

smooth=0.1

l=[]
for i in range(n):
  for c in "ACGT":
    Rex=(excounts[i].get(c,0)+exwild[i]/2+smooth/alltot[i])/(extot[i]-excounts[i].get(c,0)-exwild[i]/2+smooth)
    Rall=(allcounts[i].get(c,0)+allwild[i]/2+smooth)/(alltot[i]-allcounts[i].get(c,0)-allwild[i]/2+smooth)
    LR=log(Rex/Rall)
    l.append((LR,i,c))

l.sort(reverse=True)
for (ind,(LR,i,c)) in enumerate(l):
  if LR<0 or ind==20: break
  (aa,loc,i0)=getaa(i)
  print("%5d %c %10.1f   %s:%d"%(i+1,c,LR,aa,loc+1))

for e in range(len(header)):
  t=0
  for (LR,i,c) in l[:ind]:
    if example[e][i]==c: t+=LR
    elif example[e][i] in "ACGT": t-=LR
  print("%10.f %s"%(t,header[e]))
  
