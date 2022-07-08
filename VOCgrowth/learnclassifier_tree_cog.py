from __future__ import print_function,division
#!/usr/bin/pypy3
import sys,os
from stuff import *
from math import log

mindate="2022-05-01"
mincount=50
allowsyn=False# Allow synonymous mutations
maxleaves=20

lineages=None;numtoplin=None
if len(sys.argv)>1: mindate=sys.argv[1]
if len(sys.argv)>2:
  try:
    numtoplin=int(sys.argv[2])
  except:
    lineages=sys.argv[2].split(',')
if len(sys.argv)>3: mincount=int(sys.argv[3])
if len(sys.argv)>4: maxleaves=int(sys.argv[4])

infile='cog_metadata_sorted.csv';inputsorted=True

print('# Using input file',infile)
print("# From:",mindate)
print("# Maxleaves:",maxleaves)

allm={}
ml=[]
numl={}
for (date,lineage,mutations) in csvrows(infile,['sample_date','lineage',"mutations"]):
  if len(date)<10: continue
  if date<mindate:
    if inputsorted: break
    continue
  if lineage=="" or lineage=="Unassigned": continue
  for m in mutations.split('|'): allm[m]=allm.get(m,0)+1
  ml.append((lineage,mutations))
  numl[lineage]=numl.get(lineage,0)+1
  
print("Found",len(ml),"relevant entries since",mindate)
okm=set(m for m in allm if m!="" and allm[m]>=mincount and (allowsyn or m[:3]!="syn"))
print("Found",len(allm),"mutations, of which",len(okm),"are non-synonymous and "*(1-allowsyn)+"have occurred at least",mincount,"times")
if lineages==None: lineages=list(set(lineage for (lineage,mutations) in ml))
lineages.sort(key=lambda l: -numl[l])
if numtoplin!=None: lineages=lineages[:numtoplin]
mml=[]
for (lineage,mutations) in ml:
  if lineage in lineages: i=lineages.index(lineage)
  else: i=len(lineages)
  mml.append((i,set(mutations.split('|')).intersection(okm)))
lineages.append("Others")
print("Made reduced list")

def getstats(indexlist):
  nl=len(lineages)
  count=[0]*nl
  for i in indexlist:
    count[mml[i][0]]+=1
  s=float(sum(count))
  ent=0
  for x in count:
    if x>0: ent+=x*log(x/s)
  return (count,ent)

def splitindexlist(indexlist,mutation):
  withm=[];withoutm=[]
  for i in indexlist:
    if mutation in mml[i][1]: withm.append(i)
    else: withoutm.append(i)
  return (withm,withoutm)

class tree:
  mutation=None
  left=None# With mutation
  right=None# Without mutation
  def __init__(self,indexlist=range(len(mml)),parent=None):
    self.parent=parent
    self.indexlist=list(indexlist)
    self.count,self.ent=getstats(self.indexlist)
  def pr(self,level=0,label="Top"):
    step=4
    maxcol=16
    nl=len(lineages)
    if label=="Top":
      for j in range(min(nl,maxcol)): print(" %9s"%lineages[j],end="")
      if nl>maxcol: print(" ...",end="")
      print()
    for j in range(min(nl,maxcol)):
      print(" %9d"%self.count[j],end="")
    if nl>maxcol: print(" ...",end="")
    print(" %12.1f"%self.ent,end="  ")
    for i in range(level): print("."+" "*(step-1),end="")
    print(label)
    if self.mutation!=None:
      self.left.pr(level+1,"+"+self.mutation)
      self.right.pr(level+1,"-"+self.mutation)
  def pr2(self,mlist=[]):
    maxcol=16
    nl=len(lineages)
    if mlist==[]:
      for j in range(min(nl,maxcol)): print(" %9s"%lineages[j],end="")
      if nl>maxcol: print(" ...",end="")
      print()
    if self.mutation!=None:
      self.left.pr2(mlist+["+"+self.mutation])
      self.right.pr2(mlist+["-"+self.mutation])
      return
    for j in range(min(nl,maxcol)):
      print(" %9d"%self.count[j],end="")
    if nl>maxcol: print(" ...",end="")
    print(" %12.1f"%self.ent,end="  ")
    for m in mlist: print(m.replace("Spike","S"),end=" ")
    print()
  def split(self,mutation):
    if self.mutation!=None: raise RuntimeError("Can't split already-split node")
    (left,right)=splitindexlist(self.indexlist,mutation)
    self.mutation=mutation
    self.left=tree(left,self)
    self.right=tree(right,self)
  def getleaves(self):
    if self.mutation==None: yield self;return
    for leaf in self.left.getleaves(): yield leaf
    for leaf in self.right.getleaves(): yield leaf
  def merge(self,leaf):
    if self.mutation==None:
      tr=tree(self.indexlist+leaf.indexlist,self.parent)
      return tr,tr.ent
    else:
      l,ent_l=self.left.merge(leaf)
      r,ent_r=self.right.merge(leaf)
      p=tree(l.indexlist+r.indexlist,self.parent)
      p.mutation=self.mutation
      p.left=l
      p.right=r
      l.parent=p
      r.parent=p
      return p,ent_l+ent_r

tr=tree()
print("Total ent %.1f"%sum(leaf.ent for leaf in tr.getleaves()))
tr.pr2();print()
leaves=1
while leaves<maxleaves:
  worst=None
  for leaf in tr.getleaves():
    if worst==None or leaf.ent<worst.ent: worst=leaf
  if worst==None: break
  best=(0,None)
  for m in okm:
    (withm,withoutm)=splitindexlist(worst.indexlist,m)
    improvement=getstats(withm)[1]+getstats(withoutm)[1]-worst.ent
    #print("XXX",m,improvement)
    if improvement>best[0]: best=(improvement,m)
  if best[1]==None: print("Couldn't improve worst node");break
  worst.split(best[1])
  leaves+=1
  print("Total ent %.1f"%sum(leaf.ent for leaf in tr.getleaves()))
  tr.pr2();print();sys.stdout.flush()
print()

print("# Pruning")
while leaves>1:
  best=(-1e10,)
  for leaf in tr.getleaves():
    par=leaf.parent
    assert par!=None
    if par.left is leaf: 
      go=par.left
      keep=par.right
    else:
      go=par.right
      keep=par.left
    ent0=sum(leaf.ent for leaf in par.getleaves())
    mer,ent=keep.merge(go)
    entchg=ent-ent0
    if ent>best[0]: best=(ent,par,mer)
    #print(ent,ent0,ent-ent0)
  par.left=mer.left
  par.right=mer.right
  par.mutation=mer.mutation
  par.count=mer.count
  par.indexlist=mer.indexlist
  leaves-=1
  print("Total ent %.1f"%sum(leaf.ent for leaf in tr.getleaves()))
  tr.pr2();print();sys.stdout.flush()
