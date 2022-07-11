import sys,pickle,os
from stuff import *
from scipy.stats import norm, multivariate_normal as mvn
import numpy as np
from math import sqrt,floor,log,exp
from aliases import aliases

np.set_printoptions(precision=6,suppress=True,linewidth=200)

cachedir="cogukcachedir"
datafile='cog_metadata.csv'
Vnames=["BA.1","BA.1.1*","BA.2*"]# A name ending in '*' is considered as a prefix/ancestor, so BA.1* includes BA.1 and BA.1.17 though not BA.12
mindate=Date('2000-01-01')
maxdate=Date('2099-12-31')
mincount=5
conf=0.95
plotdots=True
plotbands=False

if len(sys.argv)>1: Vnames=sys.argv[1].split(',')
if len(sys.argv)>2: mindate=Date(sys.argv[2])
if len(sys.argv)>3: maxdate=Date(sys.argv[3])

# Valid lab locations here are UK, England, Northern_Ireland, Scotland and Wales
# NB lab location isn't necessarily the same as sample location
location="UK"

print("Labs:",location)
print("Variants considered:",' '.join(Vnames))
print("Initial date range:",mindate,"-",maxdate)
zconf=norm.ppf((1+conf)/2)
numv=len(Vnames)
cogdate=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile+'.gz')).strftime('%Y-%m-%d')

def treeclassify(mutations):
  if '|synSNP:C14599T|' in mutations and '|synSNP:C3241T|' in mutations: return "XE"
  if '|S:F486V|' in mutations:
    if '|N:P151S|' in mutations:
      if '|S:V3G|' in mutations:
        if '|S:I670V|' in mutations: return "BA.4.1.1"
        else: return "BA.4.1"
      else: return "BA.4"
    else:
      if '|synSNP:A28330G|' in mutations:
        if '|orf1ab:T5451N|' in mutations:
          return "BA.5.2"
        else:
          if '|orf1ab:V7086F|' in mutations: return "BF.1"
          else: return "BA.5.2.1"
      else:
        if '|ORF10:L37F|' in mutations: return "BA.5.1"
        else:
          if '|orf1ab:M5557I|' in mutations: return "BE.1"
          else:
            if '|S:T76I|' in mutations: return "BA.5.5"
            else:
              if '|N:E136D|' in mutations: return "BA.5.3.1"
              else:
                if '|orf1ab:R119H|' in mutations: return "BA.5.3.2"
                else: return "BA.5"
  else:
    if '|S:L452Q|' in mutations: return "BA.2.12.1"
    elif '|orf1ab:N4060S|' in mutations: return "BA.2.75"
    else:
      # Identifying "pure" BA.2 is messy due to a proliferation of BA.2.*. This rule isn't perfect, but it won't much matter as the BA.2 bit is right and it will likely get extended by COG-UK classification.
      if '|orf1ab:S135R|' in mutations and '|ORF3a:H78Y|' not in mutations and '|S:K417T|' not in mutations and '|ORF3a:L140F|' not in mutations and '|S:I68T|' not in mutations and '|orf1ab:T4175I|' not in mutations and '|S:S704L|' not in mutations and '|ORF3a:A31T|' not in mutations and '|orf1ab:S5360P|' not in mutations and '|S:F186S|' not in mutations: return "BA.2"
  return "Unassigned"

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

fn=os.path.join(cachedir,location+'_'+'_'.join(Vnames)+'__%s'%cogdate)
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    data=pickle.load(fp)
else:
  Vnames_e=[expandlin(lin) for lin in Vnames]
  data={}
  for (name,date,p2,lin,mutations) in csvrows(datafile,['sequence_name','sample_date','is_pillar_2','lineage','mutations']):
    #if p2!='Y': continue
    if location!="UK":
      country=name.split('/')[0]
      if country!=location: continue
    if not (len(date)==10 and date[:2]=="20" and date[4]=="-" and date[7]=="-"): continue
    mutations='|'+mutations+'|'
    lin_e=expandlin(lin)

    # If the COG-UK lineage is Unassigned or a prefix of tree-rule lineage, then replace it with tree-rule lineage
    if date>="2022-01-01":
      mylin=treeclassify(mutations)
      mylin_e=expandlin(mylin)
      if lin=="Unassigned" or lin_e==mylin_e[:len(lin_e)]: lin_e=mylin_e
      if mylin=="BA.2.75": lin_e=mylin_e# Special case pro tem, as COG-UK and GISAID wrongly classify this as BA.2.73
    
    # Try to assign sublineage to one of the given lineages. E.g., if Vnames=["BA.1*","BA.1.1*","BA.2"] then BA.1.14 is counted as BA.1* but BA.1.1.14 is counted as BA.1.1*
    longest=-1;ind=-1
    for (i,vn) in enumerate(Vnames_e):
      if lin_e==vn or (vn[-1]=='*' and (lin_e+'.')[:len(vn)]==vn[:-1]+'.' and len(vn)>longest): ind=i;longest=len(vn)
    if ind==-1: continue
    if date not in data: data[date]=[0]*numv
    data[date][ind]+=1
    #if date<mindate: break# If we don't abort early, then can store a cache file that works for any date range
  os.makedirs(cachedir,exist_ok=True)
  with open(fn,'wb') as fp:
    pickle.dump(data,fp)

mindate1=Date('2099-12-31')
maxdate1=Date('2000-01-01')
for date in data:
  if data[date][0]>=mincount and max(data[date][1:])>=mincount:
    mindate1=min(mindate1,Date(date))
    maxdate1=max(maxdate1,Date(date))
mindate=max(mindate,mindate1)
maxdate=min(maxdate,maxdate1)
print("Reduced date range:",mindate,"-",maxdate)

VV=[]
for date in Daterange(mindate,maxdate+1):
  v=data.get(str(date),[0]*numv)
  VV.append(v)
if VV==[]: print("No data points found");sys.exit(0)
VV=np.array(VV)

NN=VV+1e-30
n=len(NN)
if n<2: raise RuntimeError("Can't find enough samples")
# NN[timestep up to n][variant up to numv] = count

# Expand multinomial log probability about maximum:
# log(prod_r p_r ^ n_r) = constant - (1/2)*[ sum_r n_r.x_r^2 - (sum_r n_r.x_r)^2/(sum_r n_r) ] + O(x_^3)
# Note that this is invariant under x_r -> x_r + constant

# Let t = timestep, r = variant number
# Calculate, for each t, the quadratic form QF(x) = sum_r n_{t,r}x_{t,r}^2 - 1/(sum_r n_{t,r})*(sum_r n_{t,r}x_{t,r})^2
# I.e., find coefficients Q_{t,r,s} such that QF(x) = sum_{r,s} Q_{t,r,s}x_{t,r.x_{t,s}
Q=np.zeros([n,numv,numv])
Q[:,np.arange(numv),np.arange(numv)]=NN
Q-=NN[:,None,:]*NN[:,:,None]/NN.sum(axis=1)[:,None,None]

# We're looking for p_{t,r} = exp(a_r+t*b_r)/(constant independent of r), where b_r is the growth of variant r, so 
# substitute x_{t,r} = a_r + t*b_r - log(n_{t,r}), and sum over t.
# Amalgamate indexes corresponding to a_r and b_s as one single k index of size 2*numv.
# So (c_k) = (a_0,...,a_{numv-1},b_0,...,b_{numv-1}), and the quadratic (+linear) form is
#    sum_{k,l} M_{k,l}.c_k.c_l - 2*sum_k R_k.c_k + constant, i.e. c^t.M.c - 2*R.c, with MLE(c) = M^{-1}R
M=np.zeros([numv*2,numv*2])
M[:numv,:numv]=Q.sum(axis=0)
M[numv:,:numv]=M[:numv,numv:]=(Q*np.arange(n)[:,None,None]).sum(axis=0)
M[numv:,numv:]=(Q*(np.arange(n)**2)[:,None,None]).sum(axis=0)

# Add gauge-fixing terms to M, requiring the sum of a_r and sum of b_r both to be 0. (Other gauges are available.)
scale=sqrt(NN.sum())# Some kind of size scale - result doesn't depend on this unless we change it by many orders of magnitude.
M[:numv,:numv]+=scale*np.ones([numv,numv])
M[numv:,numv:]+=scale*np.ones([numv,numv])

# Calculate linear term, R
LN=np.log(NN)
R=np.zeros(numv*2)
R[:numv]=(Q*LN[:,:,None]).sum(axis=(0,1))
R[numv:]=(Q*np.arange(n)[:,None,None]*LN[:,:,None]).sum(axis=(0,1))

# Calculate provisional covariance matrix and MLE
C=np.linalg.inv(M)
c=C@R

# Calculate residuals (res), and overdispersion factor (mult), given n*(numv-1) degrees of freedom
# Effective covariance matrix, corrected for overdispersion, is then mult*C
res=c[None,:numv]+np.arange(n)[:,None]*c[None,numv:]-LN
mult=(Q*res[:,:,None]*res[:,None,:]).sum()/(n*(numv-1))
print("Residual multiplier = %.3f"%mult)

# Sampling is most convenient way of getting CrIs for crossover points
numsamp=100000
test=mvn.rvs(mean=c,cov=C*mult,size=numsamp)

# Move to gauge where baseline offset and growth are 0
test[:,:numv]-=test[:,0][:,None]
test[:,numv:]-=test[:,numv][:,None]
c[:numv]=c[:numv]-c[0]
c[numv:]=c[numv:]-c[numv]
C[:numv,:numv]=C[:numv,:numv]-C[0,:numv][None,:]-C[:numv,0][:,None]+C[0,0][None,None]
C[:numv,numv:]=C[:numv,numv:]-C[0,numv:][None,:]-C[:numv,numv][:,None]+C[0,numv][None,None]
C[numv:,:numv]=C[numv:,:numv]-C[numv,:numv][None,:]-C[numv:,0][:,None]+C[numv,0][None,None]
C[numv:,numv:]=C[numv:,numv:]-C[numv,numv:][None,:]-C[numv:,numv][:,None]+C[numv,numv][None,None]

out=[None]
for i in range(1,numv):
  print()
  grad=c[numv+i]
  graderr=sqrt(mult*C[numv+i,numv+i])*zconf
  yoff=c[i]
  t=sorted(-test[:,i]/test[:,numv+i])
  cross=t[int(numsamp/2)]
  crosslow,crosshigh=t[int(numsamp*(1-conf)/2)],t[int(numsamp*(1+conf)/2)]
  crosserr=(crosshigh-crosslow)/2
  growthstr=f"Relative growth in {Vnames[i]} vs {Vnames[0]} of {grad:.3f} ({grad-graderr:.3f} - {grad+graderr:.3f}) per day"
  doubstr=f"Doubling of ratio {Vnames[i]}/{Vnames[0]} every {log(2)/grad:.1f} ({log(2)/(grad+graderr):.1f} - {log(2)/(grad-graderr):.1f}) days"
  print(growthstr)
  print(doubstr)
  cr1=int(floor(cross))
  crossstr=f"Est'd %s/%s crossover on %s.%d +/- %.1f days"%(Vnames[i],Vnames[0],mindate+cr1,int((cross-cr1)*10),crosserr)
  print(crossstr)
  out.append((grad,graderr,yoff,cross,crosserr,growthstr,doubstr,crossstr))

datafn=location+'_%s'%('_'.join(Vnames))
visthr=2;ymin=ymax=0;ymax=-50
with open(datafn,'w') as fp:
  for date in Daterange(mindate,maxdate+1):
    t=date-mindate
    v=VV[t]
    print(date,' '.join("%6d"%x for x in v),end='',file=fp)
    mu=c[:numv]+t*c[numv:]
    var=mult*np.array([C[i,i]+2*t*C[i,numv+i]+t**2*C[numv+i,numv+i] for i in range(numv)])
    q0=mu-zconf*np.sqrt(var)
    q1=mu+zconf*np.sqrt(var)
    for i in range(numv):
      a,b=v[0]+1e-30,v[i]+1e-30
      y=log(b/a)
      prec=sqrt(a*b/(a+b))
      print(" %12g %12g"%(y,prec),end='',file=fp)
      if prec>visthr: ymin=min(ymin,y);ymax=max(ymax,y)
      #if v[0]>0 and v[i]>0: ymin=min(ymin,y);ymax=max(ymax,y)
      print(" %12g %12g"%(q0[i],q1[i]),end='',file=fp)
    print(file=fp)
print("Written data to",datafn)
  
graphfn=datafn+'.png'
ndates=maxdate-mindate+1
allothers=', '.join(Vnames[1:])
oneother=' or '.join(Vnames[1:])
number="they are" if numv>2 else "it is"

cmd=f"""
set xdata time
set key left Left reverse
set key spacing 2.5
fmt="%Y-%m-%d"
set timefmt fmt
set format x fmt
set xtics "2020-01-06", 604800
set xtics rotate by 45 right offset 0.5,0
set xtics nomirror
set grid xtics ytics lc rgb "#dddddd" lt 1
set terminal pngcairo font "sans,13" size 1728,1296
set bmargin 5.5;set lmargin 13;set rmargin 13;set tmargin 7.5
set ylabel "log_2(New {oneother} per day / New {Vnames[0]} per day)"
set style fill transparent solid 0.25
set style fill noborder

set output "{graphfn}"
set title "New cases per day in the UK of {allothers} compared with {Vnames[0]}\\nNB: This is the est'd relative growth of {allothers} compared to {Vnames[0]}, not their absolute growth. It indicates how fast {number} taking over from {Vnames[0]}\\nLarger blobs indicate more certainty (more samples). Description/caveats/current graph: http://sonorouschocolate.com/covid19/index.php/UK\\\\_variant\\\\_comparison\\nSource: Sequenced cases from COG-UK {cogdate}"
min(a,b)=(a<b)?a:b
plot [:] [{(ymin-0.5)/log(2)}:{max(ymax+0.6,1.8)/log(2)}]"""
#plot [:] [{(ymin-0.5)/log(2)}:{max(ymax+0.8*numv-0.6,1.8)/log(2)}]"""

for i in range(1,numv):
  (grad,graderr,yoff,cross,crosserr,growthstr,doubstr,crossstr)=out[i]
  if plotdots: cmd+=f""" "{datafn}" u 1:((${numv+2+4*i})/log(2)):(min(${numv+2+4*i+1},20)/{ndates/8.}) pt 5 lc {i} ps variable title "","""
  if plotbands: cmd+=f""" "{datafn}" u 1:((${numv+2+4*i+2})/log(2)):((${numv+2+4*i+3})/log(2)) w filledcurves lc {i} title "","""
  cmd+=f""" ((x/86400-{int(mindate)})*{grad}+{yoff})/log(2) lc {i} lw 2 w lines title "{growthstr}\\n{crossstr}", """
cmd+=f""" 0 lc 8 lw 3 title "{Vnames[0]} baseline" """

po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE)
p=po.stdin
p.write(cmd.encode('utf-8'))
p.close()
po.wait()
print()
print("Written graph to",graphfn)
print()

last=14
proj=28
NR=NN[-last:,:].sum(axis=0)
PR=NR/NR.sum()
stats=[]
for i in range(numv):
  grad=c[numv+i]
  graderr=sqrt(mult*C[numv+i,numv+i])*zconf
  stats.append([i,grad,graderr,PR[i],PR[i]*exp(grad*proj)])
stats=np.array(stats)
stats[:,4]=stats[:,4]/stats[:,4].sum()
orders=stats.argsort(axis=0)
types=[("original given order",0),("growth rate relative to %s"%Vnames[0],1),("relative prevalence over last %d days"%last,3),("relative prevalence projected forward %d days"%proj,4)]
for desc,col in types:
  print("Ordered by",desc)
  print("        Lineage ---------Growth---------  Last%02ddays Proj%02ddays"%(last,proj))
  for i in range(numv):
    row=stats[orders[i,col],:]
    print("%15s %6.3f (%6.3f - %6.3f)      %5.1f%%     %5.1f%%"%(Vnames[int(row[0])],row[1],row[1]-row[2],row[1]+row[2],row[3]*100,row[4]*100))
  print()
  
