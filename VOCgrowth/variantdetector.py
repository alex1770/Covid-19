import sys,time,os,pickle,argparse
from stuff import *
import numpy as np
from math import log,exp,sqrt,floor
from classify import classify, contractlin, expandlin
from scipy.optimize import minimize

np.set_printoptions(precision=6,suppress=True,linewidth=200)

mindate0=Date('2022-08-01')# Hard-coded minday for cache purposes
minmutcount=20# Ignore mutations that have occurred less than this often
cachedir='seqdatacachedir'

parser=argparse.ArgumentParser()
parser.add_argument('-f', '--mindate',        default=mindate0,     help="Min sample date of sequence")
parser.add_argument('-t', '--maxdate',        default="9999-12-31", help="Max sample date of sequence")
parser.add_argument('-a', '--effectfrom',     type=int,default=-14, help="Growth effect start date, as offset from input data datestamp")
parser.add_argument('-b', '--effectto',       type=int,default=35,  help="Growth effect end date, as offset from input data datestamp")
parser.add_argument('-g', '--gisaid',         action="store_true",  help="Use GISAID data instead of COG-UK data")
parser.add_argument('-M', '--mode',           type=int,default=0,   help="Mode: 0=evaluate individual mutations, 1=find subsets of mutations")
parser.add_argument('-m', '--givenmutations', default="",           help="Condition on this set of mutations (use a comma-separated list)")
parser.add_argument('-l', '--lineage',                              help="Condition on this lineage")
parser.add_argument('-L', '--location',       default="",           help="Location prefix")
parser.add_argument('-n', '--numdisp',        default=8,            help="Number of mutations to display horizontally")
parser.add_argument('-s', '--genomesubset',   type=int, default=2,  help="0 = Only consider mutations in a particular hand-picked subset of RBD, 1 = Only consider mutations in receptor-binding domain, 2 = Only consider mutations in spike gene, 3 = Consider mutations in any gene, but not (COG-designated) synSNPs, 4 = Consider synSNPs that are non-synonymous in some overlapping and functional ORF (e.g., A28330G), 5 = Consider any mutation, including all synSNPs")
args=parser.parse_args()

if args.mindate<mindate0: raise RuntimeError(f"Specfied mindate {args.mindate} is less that hard-coded mindate {mindate0}")

if args.gisaid:
  source='GISAID';datafile='metadata_sorted.tsv';inputsorted=True
else:
  source='COG';datafile='cog_metadata_sorted.csv';inputsorted=True

# List of overlapping ORFs for which there is evidence that they encode
# https://www.sciencedirect.com/science/article/pii/S0042682221000532
# https://virological.org/t/sars-cov-2-dont-ignore-non-canonical-genes/740/2
# 
accessorygenes=set(range(21744,21861)).union(range(25457,25580)).union(range(28284,28575))

# Hand-picked RBD subset from https://twitter.com/CorneliusRoemer/status/1576903120608600064, https://cov-spectrum.org/collections/54?region=Europe
# plus 144, 252, 484 that I added
handpickedsubset=[144,252,346,356,444,445,446,450,452,460,484,486,490,493,494]

def extractint(s):
  i=0
  while i<len(s):
    if s[i].isdigit(): break
    i+=1
  if i==len(s): return -1
  j=i
  while j<len(s):
    if not s[j].isdigit(): break
    j+=1
  return int(s[i:j])

# Determine whether mutation meets restrictions specified by -s argument
def okmut(m):
  if args.genomesubset>=5: return True
  if m[:6]=="synSNP":# Implies COG-UK
    if args.genomesubset<=3: return False
    loc=extractint(m)
    return loc in accessorygenes
  if args.genomesubset>=3: return True
  if args.gisaid and m[:6]=="Spike_": m="S:"+m[6:]
  if m[:2]!="S:": return False
  if args.genomesubset==2: return True
  loc=extractint(m)
  if args.genomesubset==1: return loc>=329 and loc<=521# RBD
  return loc in handpickedsubset

tim0=time.process_time()
datamtime=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile)).strftime('%Y-%m-%d-%H-%M-%S')
now=Date(datamtime[:10])
print("Now =",now)
print("Calculating growth from",args.mindate,"to",args.maxdate)
print("Calculating growth effect from now%+d to now%+d ="%(args.effectfrom,args.effectto),now+args.effectfrom,"to",now+args.effectto)

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
    for (dt,loc,lin,var,mutstring) in csvrows_it(fp,keys,sep=sep):
      if dt<mindate0:
        if inputsorted: break
        continue
      if len(mutstring)<=2: continue# Skip if no mutations listed
      day=datetoday(dt)
      if args.gisaid: muts=mutstring[1:-1].split(',')
      else: muts=mutstring.split('|')
      for mut in muts:
        if mut=="": print(dt,loc,lin,var,">"+mutstring+"<");oiuoiuiou
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
      linelist.append([day,loc,expandlin(classify(muts,lin,gisaid=args.gisaid)),var,l])
      for m in l: mutdaycounts[m][day]=mutdaycounts[m].get(day,0)+1
    if not inputsorted: linelist.sort(reverse=True)# Sort into reverse date order
  os.makedirs(cachedir,exist_ok=True)
  with open(fn,'wb') as fp:
    pickle.dump([linelist,num2name,name2num,mutcounts,daycounts,mutdaycounts],fp)

print("Got linelist at time",time.process_time()-tim0)

okmut_num=[okmut(x) for x in num2name]

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
  yesmuts_num={name2num[x] for x in givenmuts if x[0]!='~'}
  notmuts_num={name2num[x[1:]] for x in givenmuts if x[0]=='~'}
  if mindate!=None: minday=int(Date(mindate))
  else: minday=0
  if maxdate!=None: maxday=int(Date(maxdate))
  else: maxday=1000000
  if lineage!=None: lineage=expandlin(lineage)
  for (day,loc,lin,var,muts) in linelist:
    if day>=minday and day<maxday and yesmuts_num.issubset(muts) and notmuts_num.intersection(muts)==set() and (lineage==None or linin(lin,lineage)) and (notlineage==None or not linin(lin,notlineage)) and (location==None or loc[:len(location)]==location):
      daycounts[day]=daycounts.get(day,0)+1
      lincounts[lin]=lincounts.get(lin,0)+1
      for m in muts:
        if okmut_num[m]:
          mutdaycounts[m][day]=mutdaycounts[m].get(day,0)+1
          mutlincounts[m][lin]=mutlincounts[m].get(lin,0)+1
  return daycounts,mutdaycounts,lincounts,mutlincounts

def getgrowth(daycounts,mutdaycount):
  # log(varcount/(backgroundcount-varcount)) ~ c0+c1*(day-now) = growth*(day-crossoverday)
  if len(mutdaycount)<=3: return None
  day0=min(mutdaycount)
  day1=max(mutdaycount)
  V0=np.zeros(day1-day0+1)
  V1=np.zeros(day1-day0+1)
  for day in mutdaycount:
    x=day-mindate0
    v1=mutdaycount[day]
    v0=daycounts.get(day,0)-v1
    V0[day-day0]=v0
    V1[day-day0]=v1

  tv=V0.sum(),V1.sum()
  V0+=1e-30
  V1+=1e-30
  W=V0*V1/(V0+V1)
  if W.sum()<5: return None
  X=np.arange(day1-day0+1)
  Y=np.log(V1/V0)
  M=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
  r=np.array([sum(W*Y),sum(W*X*Y)])
  C=np.linalg.inv(M)
  c=C@r
  rho=np.exp(c[0]+c[1]*X)
  T=V0+V1;T.sort()
  # Crudely take off two degrees of freedom because rho is tuned to V0, V1 using slope and intercept. (Semi-guess, with some empirical backup.)
  mult=((V1-V0*rho)**2/rho).sum()/T[:-2].sum()
  C*=max(mult,1)
  c[0]+=c[1]*(now-day0)

  return c[0],c[1],sqrt(C[1,1]),tv

daycounts,mutdaycounts,lincounts,mutlincounts = getmutday(linelist, mindate=args.mindate, maxdate=args.maxdate, givenmuts=args.givenmutations, lineage=args.lineage, location=args.location)

print("Got all mutation-day counts at time",time.process_time()-tim0)

if args.mode==0:
  growth={}
  nsd=3
  for mut in range(nmut):
    gre=getgrowth(daycounts,mutdaycounts[mut])
    if gre!=None:
      c0,g,dg,tv=gre
      l=[]
      for gg in [g, g-nsd*dg, g+nsd*dg]:
        if g*gg<0: gg=0
        R0=exp(c0+gg*args.effectfrom)
        R1=exp(c0+gg*args.effectto)
        greffect=gg*(R1-R0)/((1+R0)*(1+R1))
        l.append(greffect)
      if l[1]>l[2]: l=[l[0],l[2],l[1]]
      growth[mut]=g,dg,l[0],l[1],l[2],tv
  l=list(growth)
  #l.sort(key=lambda x:-(growth[x][0]-nsd*growth[x][1]))
  #l.sort(key=lambda x:-abs((growth[x][0]-gr0)/growth[x][1]))
  l.sort(key=lambda x:-growth[x][3])

  print("     Mutation          -------- %Growth ---------     ---- %Growth effect ----- Gr-signif NumNonVar  NumVar  Examples")
  nm=0
  for mut in l:
    g,dg,greff,greff_low,greff_high,tv=growth[mut]
    (gl,gh)=(g-nsd*dg,g+nsd*dg)
    if gl<=0 and gh>=0: continue
    
    # Growth
    print("%-20s  %7.2f (%7.2f - %7.2f)"%(num2name[mut],g*100,gl*100,gh*100),end='')
    
    # Growth effect
    print("   %7.2f (%7.2f - %7.2f)"%(greff*100,greff_low*100,greff_high*100),end='')
    
    print("   %7.2f   %7d %7d"%(greff/((greff_high-greff_low)/(2*nsd)),tv[0],tv[1]),end='')
    
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
      elif n0==0 or n/n0<0.05: break
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



if args.mode==1:
  
  mlist=[m for m in range(nmut) if okmut_num[m]]
  mlist.sort(key=lambda m:(extractint(num2name[m]),num2name[m]))
  mindate=Date(args.mindate)
  maxdate=min(args.maxdate,Date(linelist[0][0]))
  ndays=maxdate-mindate+1
  
  linelist1=[]
  for x in linelist:
    if x[0]>maxdate: continue
    if x[0]<mindate: break
    linelist1.append(x)
    x[4]=[m for m in x[4] if okmut_num[m]]
  linelist=linelist1
  print("Filtered linelist at time",time.process_time()-tim0)
  print()
  
  def LL(xx):
    gg=bit@xx[1<<n:]
    LL1=nn0@xx[:1<<n]+nn1@gg
    den=np.exp(xx[:1<<n,None]+gg[:,None]*np.arange(ndays)).sum(axis=0)
    LL1-=nn@np.log(den)
    return LL1

  def dLL(xx):
    dLL1=np.concatenate([nn0,nn1@bit])
    gg=bit@xx[1<<n:]
    dens=np.exp(xx[:1<<n,None]+gg[:,None]*np.arange(ndays))
    den=dens.sum(axis=0)
    dLL1[:1<<n]-=((dens*nn)/den).sum(axis=1)
    dLL1[1<<n:]-=((bit.T@dens)*nnt/den).sum(axis=1)
    return dLL1

  def NLL(xx): return -LL(xx)/(((1<<n)+n)*ndays)
  def NdLL(xx): return -dLL(xx)/(((1<<n)+n)*ndays)

  def LLpr(xx):
    gg=bit@xx[1<<n:]
    dens=np.exp(xx[:1<<n,None]+gg[:,None]*np.arange(ndays))
    den=dens.sum(axis=0)
    print("                                   ",end="")
    for I in range(1<<n): print(" %10.5f"%xx[I],end="")
    print()
    print("                                   ",end="")
    for I in range(1<<n): print(" %10.5f"%gg[I],end="")
    print()
    lltot=0
    for t in range(ndays):
      print("%3d |"%t,end="")
      for I in range(1<<n): print(" %6d"%nnx[I,t],end="")
      print(" |",end="")
      for I in range(1<<n): print(" %10.4g"%(dens[I][t]/den[t]),end="")
      ll=0
      for I in range(1<<n): ll+=nnx[I,t]*log(dens[I][t]/den[t])
      print(" |",end="")
      for I in range(1<<n): print(" %7.5f"%(nnx[I,t]/nn[t]),end="")
      lltot+=ll
      print(" | %14.6f %14.6f"%(ll,lltot))
    assert abs(LL(xx)-lltot)<1e-3

  # Growth effect, defined up to an additive constant
  # GE(xx,t1)-GE(xx,t0) is well-defined
  def GE(xx,t):
    gg=bit@xx[1<<n:]
    den=np.exp(xx[:1<<n]+gg*t)
    return gg@den/den.sum()
    
  M=[]
  thr=0.005# Require at least this much improvement in growth effect
  ge0=0
  best0=[0,]
  while 1:
    Mints=[extractint(num2name[m]) for m in M]
    deb=0
    best=(-1,)
    for mnew in mlist:
      if extractint(num2name[mnew]) in Mints: continue
      # Contemplating M -> M u {m}
      
      M1=M+[mnew]
      n=len(M1)
      # a_I, g_i
      # I <-> 2^n
      # Use gauge: a_{empty}=0; (g_i don't need gauge fixing)
      
      # n_{t,I} = number of instances of mutation pattern I on day t
      # nn0[I] = sum_t n_{t,I}
      # nn1[I] = sum_t t*n_{t,I}
      # nn[t]  = sum_I n_{t,I}
      nn0=np.zeros(1<<n)
      nn1=np.zeros(1<<n)
      nn=np.zeros(ndays)
      nnt=np.zeros(ndays)
      if deb: nnx=np.zeros([1<<n,ndays])
      for x in linelist:
        I=0
        for ml in x[4]:
          if ml in M1: I+=1<<M1.index(ml)
        t=x[0]-mindate
        nn0[I]+=1
        nn1[I]+=t
        nn[t]+=1
        nnt[t]+=t
        if deb: nnx[I,t]+=1
  
      bit=np.zeros([1<<n,n])
      for i in range(1<<n):
        for j in range(n):
          bit[i,j]=(i>>j)&1
  
      bounds=[(-20,10)]*(1<<n)+[(-0.5,0.5)]*n
      bounds[0]=(0,0)
      for i in range(1,1<<n):
        if nn0[i]==0: bounds[i]=(-30,-30)
      xx=[0]*((1<<n)+n)
      res=minimize(NLL,xx,bounds=bounds, jac=NdLL, method="SLSQP", options={'ftol':1e-20, 'maxiter':10000})
      xx=res.x
      if not res.success: raise RuntimeError(res.message)
      err=0
      for i in range(len(xx)):
        if bounds[i][0]<bounds[i][1] and (xx[i]<bounds[i][0]+1e-3 or xx[i]>bounds[i][1]-1e-3):
          err=1
          if i<(1<<n): print("Error: intercept of",'+'.join(num2name[M1[j]] for j in range(n) if bit[i,j]),"hit bound")
          else: print("Error: growth of",num2name[M1[i-(1<<n)]],"hit bound")
      if err: raise RuntimeError("Optimisation hit bounds")
      for m in M1:
        print("%15s"%num2name[m],end="")
      ge=GE(xx,now+args.effectto-mindate)-GE(xx,now+args.effectfrom-mindate)
      #print(" ",xx,nn0,"| %7.5f"%ge)
      for i in range(n): print("  %7.4f"%xx[(1<<n)+i],end="")
      print(" | %7.4f"%ge)
      if ge>best[0]: best=(ge,mnew,xx)
      if deb: LLpr(xx)
    ge,m,xx=best
    print("Best growth effect",ge,"using",num2name[m])
    print()
    if ge-best0[0]<thr: break
    best0=best
    M.append(m)
  print("Final choice")
  ge,m,xx=best0
  n=len(M)
  for i,m in enumerate(M):
    print("%15s  %7.4f"%(num2name[m],xx[(1<<n)+i]))
  print("Growth effect %.4f"%ge)
