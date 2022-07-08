from __future__ import print_function,division
#!/usr/bin/pypy3
import sys,os
from stuff import *
from math import log

mindate="2022-05-01"
mincount=50
allowsyn=False# Allow synonymous mutations
maxleaves=50

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
    maxcol=30
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
    maxcol=30
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
  # Merge indexlist into the tree self, returning a new allocated tree
  def merge(self,indexlist):
    if self.mutation==None:
      tr=tree(self.indexlist+indexlist,self.parent)
      return tr
    else:
      left,right=splitindexlist(indexlist,self.mutation)
      l=self.left.merge(left)
      r=self.right.merge(right)
      p=tree(l.indexlist+r.indexlist,self.parent)
      p.mutation=self.mutation
      p.left=l
      p.right=r
      l.parent=p
      r.parent=p
      return p
  def check(self):
    count,ent=getstats(self.indexlist)
    if count!=self.count or abs(ent-self.ent)>1e-6: return False
    if self.mutation==None:
      return self.left==None and self.right==None
    if self.left==None or self.right==None: return False
    if self.left.parent!=self or self.right.parent!=self: return False
    return self.left.check() and self.right.check()
  def leafent(self):
    return sum(leaf.ent for leaf in self.getleaves())

tr=tree()
print("Total ent %.1f"%tr.leafent())
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
  print("Total ent %.1f"%tr.leafent())
  tr.pr2();print();sys.stdout.flush()
print()

print("# Pruning")
print()
while leaves>1:
  best=(-1e10,)
  if not tr.check(): raise RuntimeError("A")
  for leaf in tr.getleaves():
    par=leaf.parent
    assert par!=None
    if par.left is leaf: 
      go=par.left
      keep=par.right
    else:
      go=par.right
      keep=par.left
    mer=keep.merge(go.indexlist)
    if not mer.check(): raise RuntimeError("B")
    entchg=mer.leafent()-par.leafent()
    if entchg>best[0]: best=(entchg,par,mer)
  entchg,par,mer=best
  par.left=mer.left
  par.right=mer.right
  par.mutation=mer.mutation
  if par.mutation is not None:
    par.left.parent=par
    par.right.parent=par
  par.count=mer.count
  par.indexlist=mer.indexlist
  if not par.check(): raise RuntimeError("C")
  if not tr.check(): raise RuntimeError("D")
  leaves-=1
  print("Total ent %.1f"%tr.leafent())
  tr.pr2();print();sys.stdout.flush()
