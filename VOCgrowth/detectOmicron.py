import sys,time,os,pickle
from stuff import *
import numpy as np
from math import log,sqrt,floor

mindate='2021-01-01'
maxdate='9999-12-31'

datafile='cog_metadata_sorted.csv'

# Δ69-70 = del_21765_6 in COG notation
numut0=['S:'+x for x in 'A67V, T95I, G142D, L212I, ins214EPE, G339D, S371L, S373P, S375F, K417N, N440K, G446S, S477N, T478K, E484A, Q493K, G496S, Q498R, N501Y, Y505H, T547K, D614G, H655Y, N679K, P681H, N764K, D796Y, N856K, Q954H, N969K, L981F'.replace(',','').split()]
numut0+='E:T9I M:Q19E M:A63T N:P13L N:R203K N:G204R'.split()
numut0+='del6970'
numut={}

mutcounts={}
n=0
with open(datafile,'r') as fp:
  for (dt,lin,var,del6970,mut) in csvrows_it(fp,['sample_date','lineage','scorpio_call','del_21765_6','mutations']):
    if dt>maxdate: continue
    if dt<mindate: break
    day=datetoday(dt)
    muts=mut.split('|')
    if del6970=='del': muts.append('del6970')
    for mut in muts:
      mutcounts[mut]=mutcounts.get(mut,0)+1
    n+=1

for mut in numut0:
  numut[mut]=-log((mutcounts.get(mut,0)+1)/(n+1))/len(numut0)

with open(datafile,'r') as fp:
  for (name,dt,lin,var,del6970,mut) in csvrows_it(fp,['sequence_name','sample_date','lineage','scorpio_call','del_21765_6','mutations']):
    if dt>maxdate: continue
    if dt<mindate: break
    day=datetoday(dt)
    muts=mut.split('|')
    if del6970=='del': muts.append('del6970')
    sc0=sc1=0
    for mut in muts:
      if mut in numut: sc0+=1;sc1+=numut[mut]
    print("%4d %9.3f %12s %20s"%(sc0,sc1,dt,name),lin)
