import sys,os
from stuff import *
from scipy.stats import binom

c='Europe / United Kingdom'
mindate='2021-01-01'

#VV=['B.1.1.7','B.1.617','B.1.617.1','B.1.617.2','B.1.617.3']
#VV=['B.1.1.7','B.1.617.2','AY.1','AY.2','AY.3','AY.4','AY.5','AY.6','C.1.2']
#VV=['B.1.617.2','AY.1','AY.2','AY.3','AY.4','AY.5','AY.6','AY.4.2']
#VV=['B.1.617.2','AY.4','AY.4.2','AY.43']
#VV=['B.1.617.2','AY.4','AY.4.2','AY.43','B.1.1.529']
#VV=['B.1','B.1.2','B.1.243','B.1.1.7','B.1.1.519','B.1.427','B.1.429']
#VV=['B.1','B.1.2','B.1.243','B.1.1.7','B.1.1.519','B.1.427','B.1.429']
#VV=['B.1.617.2','AY.*','BA.1','BA.1.1','BA.2','BA.3']
VV=['Spike_F486V', 'N_P151S', 'NSP8_N118S', 'BA.2*', 'BA.3', 'BA.4*', 'BA.5*', 'Unassigned']

if len(sys.argv)>1: c0=sys.argv[1]
if len(sys.argv)>2: c1=sys.argv[2]
if len(sys.argv)>3: mindate=sys.argv[3]

try:
  t0=os.path.getmtime('metadata.tsv')
except FileNotFoundError:
  t0=-1e30

try:
  t1=os.path.getmtime('metadata_sorted.tsv')
except FileNotFoundError:
  t1=-1e30

if t0<0 and t1<0: raise FileNotFoundError("Could not find GISAID files metadata.tsv or metadata_sorted.tsv")
if t1>=t0:
  infile='metadata_sorted.tsv';inputsorted=True
else:
  infile='metadata.tsv';inputsorted=False

print('Using input file',infile)
print("Comparing location",c0)
print("with location",c1)
print("from",mindate)

num0={};num1={}
t0=t1=0
for (date,loc,lineage,mutations) in csvrows(infile,['Collection date','Location','Pango lineage','AA Substitutions'],sep='\t'):
  in0=(loc[:len(c0)]==c0)
  in1=(loc[:len(c1)]==c1)
  if not (in0 or in1): continue
  if len(date)!=10 or date[:2]!="20": continue
  if date<mindate:
    if inputsorted: break
    continue
  mlist=mutations[1:-1].split(',')
  if in0: t0+=1
  if in1: t1+=1
  for m in mlist:
    if in0: num0[m]=num0.get(m,0)+1
    if in1: num1[m]=num1.get(m,0)+1

print("Found",t0,"sequences and",len(num0),"distinct mutations in",c0)
print("Found",t1,"sequences and",len(num1),"distinct mutations in",c1)

def CI(n,k,conf=0.99):
  eps=1e-3
  # find p0,p1 such that
  # P(B(n,p0)< k)=(1-conf)/2
  # P(B(n,p1)<=k)=(1+conf)/2
  a=0;b=1
  while 1:
    c=(a+b)/2
    if b-a<eps: break
    if binom.cdf(k,n,c)<(1-conf)/2: b=c
    else: a=c
  p1=c
  a=0;b=1
  while 1:
    c=(a+b)/2
    if b-a<eps: break
    if binom.cdf(k-1,n,c)<(1+conf)/2: b=c
    else: a=c
  p0=c
  return p0,p1

both=set(num0).union(num1)
l=[]
for m in both:
  n0=num0.get(m,0)
  n1=num1.get(m,0)
  # Is n0 out of t0 significantly different from n1 out of t1?
  # Actually, the better question is what is the minimum significant ratio change?
  pmin0,pmax0=CI(t0,n0)
  pmin1,pmax1=CI(t1,n1)
  rchg=min(pmax1/pmin0,pmax0/pmin1)
  if rchg<1: l.append((rchg,m,n0,t0,n1,t1))

l.sort()
for (pv,m,n0,t0,n1,t1) in l:
  print("%12g %16s %8d %8d %8d %8d"%(pv,m,n0,t0,n1,t1))
