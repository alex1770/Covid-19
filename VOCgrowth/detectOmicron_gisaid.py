# Do this so can use pypy
from __future__ import division, print_function

import sys,pickle
from stuff import *
from math import log

mindate='2021-09-01'
maxdate='9999-12-31'

datafile='gisaid_metadata.tsv'
thr=0.4
if len(sys.argv)>1: mindate=sys.argv[1]
if len(sys.argv)>2: datafile=sys.argv[2]
if len(sys.argv)>3: thr=float(sys.argv[3])
print("Dates",mindate,"-",maxdate)
print("Using data from",datafile)
print("Using threshold",thr)

# Example of Omicron:
# (NSP5_P132H,Spike_H69del,Spike_T95I,Spike_A67V,Spike_S373P,Spike_H655Y,N_R203K,Spike_N969K,Spike_N856K,Spike_G142D,NSP3_A1892T,Spike_Q954H,N_P13L,NSP3_L1266I,N_R32del,M_Q19E,Spike_N440K,NSP4_T492I,NSP6_L105del,Spike_N679K,Spike_N764K,Spike_L212I,NSP6_G107del,NSP6_I189V,Spike_T547K,M_D3G,Spike_D796Y,N_G204R,Spike_V143del,M_A63T,Spike_K417N,NSP6_S106del,Spike_S371L,Spike_G339D,NSP3_S1265del,NSP14_I42V,Spike_P681H,Spike_Y144del,Spike_ins214EPE,N_S33del,Spike_G446S,N_E31del,NSP3_K38R,Spike_N211del,E_T9I,Spike_V70del,Spike_L981F,NSP12_P323L,Spike_Y145del,Spike_D614G)

numut0=['Spike_'+x for x in 'A67V T95I G142D L212I ins214EPE G339D S371L S373P S375F K417N N440K G446S S477N T478K E484A Q493K G496S Q498R N501Y Y505H T547K D614G H655Y N679K P681H N764K D796Y N856K Q954H N969K L981F'.split()]
numut0+='NSP3_K38R NSP3_V1069I NSP3_L1266I NSP3_A1892T NSP4_T492I NSP5_P132H NSP6_L105del NSP6_G107del NSP6_A189V NSP12_P323L NSP14_I42V E_T9I M_D3G M_Q19E M_A63T N_P13L N_E31del N_R203K N_G204R'.split()
numut={}

mutcounts={}
n=0
with open(datafile,'r') as fp:
  for (dt,lin,mut) in csvrows_it(fp,['Collection date','Pango lineage','AA Substitutions'],sep='\t'):
    if dt>maxdate: continue
    if dt<mindate: continue
    muts=mut[1:-1].split(',')
    for mut in muts:
      mutcounts[mut]=mutcounts.get(mut,0)+1
    n+=1

for mut in numut0:
  numut[mut]=-log((mutcounts.get(mut,0)+1)/(n+1))/len(numut0)

l=[]
with open(datafile,'r') as fp:
  for (name,dt,lin,mut) in csvrows_it(fp,['Virus name','Collection date','Pango lineage','AA Substitutions'],sep='\t'):
    if dt>maxdate: continue
    if dt<mindate: continue
    muts=mut[1:-1].split(',')
    sc0=sc1=0
    for mut in muts:
      if mut in numut: sc0+=1;sc1+=numut[mut]
    if sc1>thr: l.append((sc1,sc0,dt,name,lin))

l.sort()
for (sc1,sc0,dt,name,lin) in l:
  print("%4d %9.3f %12s %20s"%(sc0,sc1,dt,name),lin)
