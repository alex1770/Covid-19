import sys,time,os,pickle,argparse
from stuff import *
import numpy as np
from math import log,sqrt,floor
from variantaliases import aliases

mindate0=Date('2022-05-01')# Hard-coded minday
minmutcount=50# Ignore mutations that have occurred less than this often
cachedir='seqdatacachedir'

parser=argparse.ArgumentParser()
parser.add_argument('-f', '--mindate',        default=mindate0,     help="Min sample date of sequence")
parser.add_argument('-t', '--maxdate',        default="9999-12-31", help="Max sample date of sequence")
parser.add_argument('-g', '--gisaid',         action="store_true",  help="Use GISAID data instead of COG-UK data")
parser.add_argument('-m', '--givenmutations', default="",           help="Condition on this set of mutations (use a comma-separated list)")
parser.add_argument('-l', '--lineage',                              help="Condition on this lineage")
parser.add_argument('-L', '--location',       default="",           help="Location prefix")
parser.add_argument('-n', '--numdisp',        default=8,            help="Number of mutations to display horizontally")
args=parser.parse_args()

if args.mindate<mindate0: raise RuntimeError(f"Specfied mindate {args.mindate} is less that hard-coded mindate {mindate0}")

if args.gisaid:
  source='GISAID';datafile='metadata_sorted.tsv';inputsorted=True
else:
  source='COG';datafile='cog_metadata_sorted.csv';inputsorted=True

ecache={}
def expandlin(lin):
  if lin in ecache: return ecache[lin]
  for (short,long) in aliases:
    s=len(short)
    if lin[:s+1]==short+".": ecache[lin]=long+lin[s:];return ecache[lin]
  ecache[lin]=lin
  return lin

ccache={}
def contractlin(lin):
  if lin in ccache: return ccache[lin]
  lin=expandlin(lin)
  for (short,long) in aliases:
    l=len(long)
    if lin[:l+1]==long+".": ccache[lin]=short+lin[l:];return ccache[lin]
  ccache[lin]=lin
  return lin

tim0=time.process_time()
datamtime=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile)).strftime('%Y-%m-%d-%H-%M-%S')
id=source+'_'+datamtime+'_'+mindate0+'_'+str(minmutcount)
fn=os.path.join(cachedir,id)
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    (linelist,num2name,name2num,mutcounts,daycounts,mutdaycounts)=pickle.load(fp)
  nmut=len(num2name)
else:
  print("Reading %s"%datafile)
  if args.gisaid: keys=["Collection date","Location","Pango lineage","Variant","AA Substitutions"];sep='\t'
  else: keys=['sample_date',"adm1",'lineage','scorpio_call','mutations'];sep=','
  mutcounts={}
  with open(datafile,'r') as fp:
    for (dt,loc,lin,var,mut) in csvrows_it(fp,keys,sep=sep):
      if dt<mindate0:
        if inputsorted: break
        continue
      day=datetoday(dt)
      if args.gisaid: muts=mut[1:-1].split(',')
      else: muts=mut.split('|')
      for mut in muts:
        mutcounts[mut]=mutcounts.get(mut,0)+1
  num2name=[]
  name2num={}
  for mut in sorted(mutcounts):
    if mutcounts[mut]>=minmutcount:
      name2num[mut]=len(num2name)
      num2name.append(mut)
  print("Found %d mutations of which %d occur at least %d times, in %.3fs"%(len(mutcounts),len(num2name),minmutcount,time.process_time()-tim0))
  nmut=len(num2name)
  linelist=[]
  mutcounts=[mutcounts[mut] for mut in num2name]
  daycounts={}
  mutdaycounts=[{} for m in range(nmut)]
  with open(datafile,'r') as fp:
    for (dt,loc,lin,var,mut) in csvrows_it(fp,keys,sep=sep):
      if dt<mindate0:
        if inputsorted: break
        continue
      day=datetoday(dt)
      daycounts[day]=daycounts.get(day,0)+1
      if args.gisaid: muts=mut[1:-1].split(',')
      else: muts=mut.split('|')
      l=[]
      for mut in muts:
        if mut in name2num: l.append(name2num[mut])
      linelist.append([day,loc,expandlin(lin),var,l])
      for m in l: mutdaycounts[m][day]=mutdaycounts[m].get(day,0)+1
    if not inputsorted: linelist.sort(reverse=True)# Sort into reverse date order
  os.makedirs(cachedir,exist_ok=True)
  with open(fn,'wb') as fp:
    pickle.dump([linelist,num2name,name2num,mutcounts,daycounts,mutdaycounts],fp)

print("Got linelist at time",time.process_time()-tim0)

# Does 'lin' fit in with the (possibly wildcarded) 'lineage'?
# Note that BA.1.2 is considered part of BA.1*, but
#           BA.12 is not considered part of BA.1*,
# so it's not just a prefix test.
# Input is assumed to be in expanded form.
def linin(lin,lineage):
  if lin==lineage: return True
  if lineage[-1]!='*': return False
  if lin==lineage[:-1]: return True
  return lin[:len(lineage)]==lineage[:-1]+'.'
  

def getmutday(linelist,mindate=None,maxdate=None,givenmuts=[],lineage=None,notlineage=None,location=None):# Might need a givennot argument too
  daycounts={}
  mutdaycounts=[{} for m in range(nmut)]
  mutlincounts=[{} for m in range(nmut)]
  lincounts={}
  if type(givenmuts)==str: givenmuts=(givenmuts.split(',') if givenmuts!="" else [])
  givenmuts_num={name2num[x] for x in givenmuts}
  if mindate!=None: minday=int(Date(mindate))
  else: minday=0
  if maxdate!=None: maxday=int(Date(maxdate))
  else: maxday=1000000
  if lineage!=None: lineage=expandlin(lineage)
  for (day,loc,lin,var,muts) in linelist:
    if day>=minday and day<maxday and givenmuts_num.issubset(muts) and (lineage==None or linin(lin,lineage)) and (notlineage==None or not linin(lin,notlineage)) and (location==None or loc[:len(location)]==location):
      daycounts[day]=daycounts.get(day,0)+1
      lincounts[lin]=lincounts.get(lin,0)+1
      for m in muts:
        mutdaycounts[m][day]=mutdaycounts[m].get(day,0)+1
        mutlincounts[m][lin]=mutlincounts[m].get(lin,0)+1
  return daycounts,mutdaycounts,lincounts,mutlincounts

def getgrowth(daycounts,mutdaycount,invert=False):
  # log(varcount/(backgroundcount-varcount)) ~ c0+c1*(day-minday) = growth*(day-crossoverday)
  m=np.zeros([2,2])
  r=np.zeros(2)
  s0=syy=0
  tv0=tv1=0
  for day in mutdaycount:
    x=day-mindate0
    v1=mutdaycount[day]
    v0=daycounts.get(day,0)-v1
    if invert: (v0,v1)=(v1,v0)
    if v0>0 and v1>0:
      y=log(v1/v0)
      w=1/(1/v0+1/v1)
      m[0,0]+=w
      m[0,1]+=w*x
      m[1,0]+=w*x
      m[1,1]+=w*x*x
      r[0]+=w*y
      r[1]+=w*x*y
      s0+=1
      syy+=w*y*y
      tv0+=v0;tv1+=v1
  if m[0,0]<5: return None#(0,1000),(0,10),(tv0,tv1)#alter
  mi=np.linalg.pinv(m)
  c=mi@r#c=np.linalg.solve(m,r)
  cv=[mi[0,0],mi[1,1]]# These should be the variances of c[0],c[1]

  # alter - not using residual yet
  # Want sum( w*(c0+c1*x-y)^2 )
  # vr = (c[0]**2*m[0,0] + 2*c[0]*c[1]*m[0,1] - 2*c[0]*r[0] + c[1]**2*m[1,1] - 2*c[1]*r[1] + syy)/s0
  #print(vr)
  # Investigate simple correction for overdispersion. Not sure it's right yet.
  # Try crossing number stats

  return (c[0],sqrt(cv[0])),(c[1],sqrt(cv[1])),(tv0,tv1)
  # This form is nicer to interpret (and minday-independent), but will become singular if c[1]=0:
  # return (mindate0-c[0]/c[1],sqrt(cv[0])/c[1]),(c[1],sqrt(cv[1]))

#daycounts,mutdaycounts,lincounts=getmutday(linelist)
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-07-01',maxdate='2021-10-01')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,givenmuts={'S:Y145H'})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-07-01',maxdate='2021-08-16')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-07-01',maxdate='2021-08-01',givenmuts={'S:G142D'})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-01-01',maxdate='2021-08-01',lineage='AY.4')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-01-01',maxdate='2021-08-15',givenmuts={'S:T95I'},lineage='AY.4')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-07-01',maxdate='2021-08-15',givenmuts={'S:T95I'})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-01-01',maxdate='2021-08-01',lineage='AY.4')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-09-01',maxdate='2021-10-07',givenmuts={'S:Y145H'})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-08-01',lineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-10-14',lineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-10-14',notlineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-08-01',maxdate='2021-10-10',lineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-10-14',notlineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-08-01',maxdate='2021-10-10',notlineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-10-18',notlineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-08-01',notlineage='AY.4.2',maxdate='2021-10-12')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-08-01',notlineage='AY.4.2',givenmuts={'N:Q9L'})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-08-01',notlineage='AY.4.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-12-25',lineage='BA.1')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2021-12-25',lineage='BA.1',givenmuts={'S:N440K'})
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2022-01-17',lineage='BA.1')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2022-01-17',lineage='BA.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2022-03-01',lineage='BA.2')
#daycounts,mutdaycounts,lincounts=getmutday(linelist,mindate='2022-03-01')
#daycounts,mutdaycounts,lincounts,mutlincounts=getmutday(linelist,mindate='2022-06-01',givenmuts={'S:F486V'})
#daycounts,mutdaycounts,lincounts,mutlincounts=getmutday(linelist,mindate='2022-06-01',givenmuts={'Spike_F486V'})
#daycounts,mutdaycounts,lincounts,mutlincounts=getmutday(linelist,mindate='2022-06-01',givenmuts={'Spike_F486V'},location="Europe / Denmark")
#daycounts,mutdaycounts,lincounts,mutlincounts=getmutday(linelist,mindate='2022-06-01',givenmuts={'Spike_F486V'},location="Europe / Austria / Vienna")
#daycounts,mutdaycounts,lincounts,mutlincounts=getmutday(linelist,mindate='2022-05-01')
#daycounts,mutdaycounts,lincounts,mutlincounts=getmutday(linelist,mindate='2022-05-01',location="Asia / India")
daycounts,mutdaycounts,lincounts,mutlincounts = getmutday(linelist, mindate=args.mindate, maxdate=args.maxdate, givenmuts=args.givenmutations, lineage=args.lineage, location=args.location)

print("Got all mutation-day counts at time",time.process_time()-tim0)

if 1:
  growth={};tv={}
  okmuts=[]
  for mut in range(nmut):
    gr=getgrowth(daycounts,mutdaycounts[mut])
    if gr!=None:
      growth[mut]=gr[1]
      tv[mut]=gr[2]
      okmuts.append(mut)

  with open('tempvargr','w') as fp:
    sg=[growth[x][0]/growth[x][1] for x in growth]
    n=len(sg)
    d=int(sqrt(n))
    sg.sort()
    x0=sg[int(0.01*n)]
    x1=sg[-max(int(0.01*n),1)]
    hist=[0]*d
    low=high=0
    for x in sg:
      i=int(floor((x-x0)/(x1-x0)*d))
      if i<0: low+=1
      elif i<d: hist[i]+=1
      else: high+=1
    print("# Low %d"%low,file=fp)
    for (i,k) in enumerate(hist):
      print("%8.3f  %8d"%(x0+(i+.5)/d*(x1-x0),hist[i]),file=fp)
    print("# High %d"%high,file=fp)

  sd=8
  gr0=0.0
  def gval(mut):
    if tv[mut][0]+tv[mut][1]==0: return -1e9
    gr=growth[mut]
    p=tv[mut][0]/(tv[mut][0]+tv[mut][1])
    g,dg=gr[0],gr[1]
    if g<0: g=-g;p=1-p
    return (g-sd*dg)*p
  #l=[mut for mut in growth if growth[mut][0]>0]
  l=list(growth)
  #l.sort(key=lambda x:-(growth[x][0]-sd*growth[x][1]))
  #l.sort(key=lambda x:-abs((growth[x][0]-gr0)/growth[x][1]))
  l.sort(key=lambda x:-gval(x))

  print("     Mutation          -------- %Growth ---------     ---- %Growth effect -----    GE-1sd Gr-signif NumNonVar  NumVar Examples")
  nm=0
  for mut in l:
    gr=growth[mut]
    (g,gl,gh)=(gr[0],gr[0]-sd*gr[1],gr[0]+sd*gr[1])
    #if gl<gr0: break
    #if abs(gr[0])<sd*gr[1]: break
    print("%-20s  %7.3f (%7.3f - %7.3f)"%(num2name[mut],g*100,gl*100,gh*100),end='')
    p=tv[mut][0]/(tv[mut][0]+tv[mut][1])
    if g>0: (ga,gla,gha)=(g*p,gl*p,gh*p)
    else: (ga,gla,gha)=(-g*(1-p),-gh*(1-p),-gl*(1-p))
    print("   %7.3f (%7.3f - %7.3f)"%(ga*100,gla*100,gha*100),end='')
    print("   %7.2f"%((ga-.5*(gha-gla))*100),end='')
    print("   %7.2f   %7d %7d"%((gr[0]-gr0)/gr[1],tv[mut][0],tv[mut][1]),end='')
    ml=list(mutlincounts[mut])
    score={}# score = number of lineages with mutation if it's a growing mutation, or number without mutation if it's a falling mutation
    for lin in ml:
      if g>0: score[lin]=mutlincounts[mut][lin]
      else: score[lin]=lincounts[lin]-mutlincounts[mut][lin]
    ml.sort(key=lambda lin:-score[lin])
    n0=None
    for lin in ml[:3]:
      n=score[lin]
      if n0==None: n0=n
      elif n/n0<0.05: break
      print("  %s:%d"%(contractlin(lin),n),end="")
    print()
    nm+=1
    if nm==30: break
  
  nmd=min(nm,args.numdisp)
  print()
  days=set()
  for mut in l[:nmd]: days.update(mutdaycounts[mut])
  print("                All",end='')
  for mut in l[:nmd]: print(" %20s"%num2name[mut],end='')
  print()
  days=sorted(days)
  for day in days:
    v0=daycounts[day]
    print(daytodate(day),"  %6d"%v0,end='')
    for mut in l[:nmd]:
      vm=mutdaycounts[mut].get(day,0)
      v1=v0-vm
      if vm>0 and v1>0: print("         %6.3f"%(log(vm/v1)),end='')
      else: print("              -",end='')
      print(" %5d"%vm,end='')
    print()

