import sys,pickle,os,argparse
from stuff import *
from scipy.stats import norm, multivariate_normal as mvn
from scipy.optimize import minimize
from scipy.special import gammaln
import numpy as np
from math import sqrt,floor,log,exp
from variantaliases import aliases

np.set_printoptions(precision=6,suppress=True,linewidth=200)

cachedir="cogukcachedir"
datafile='cog_metadata.csv'
conf=0.95

parser=argparse.ArgumentParser()
parser.add_argument('-c', '--mincount',    type=int,default=5,    help="Minimum variant count considered")
parser.add_argument('-d', '--DM',          action="store_true",   help="Use Dirichlet-Multinomial regression instead of logistic approximation")
parser.add_argument('-f', '--mindate',     default="2022-01-01",  help="Min sample date of sequence")
parser.add_argument('-t', '--maxdate',     default="9999-12-31",  help="Max sample date of sequence")
parser.add_argument('-l', '--lineages',    default="BA.4*,BA.5*", help="Comma-separated list of lineages/variants")
parser.add_argument('-p', '--plotpoints',  action="store_true",   help="Whether to plot points corresponding to log(num(variant)/num(base variant))")
parser.add_argument('-b', '--plotbands',   action="store_true",   help="Whether to plot confidence bands around best-fit lines")
args=parser.parse_args()

Vnames=args.lineages.split(',')
mindate=Date(args.mindate)
maxdate=Date(args.maxdate)

# Valid lab locations here are UK, England, Northern_Ireland, Scotland and Wales
# NB lab location isn't necessarily the same as sample location
location="UK"

print("Labs:",location)
print("Variants considered:",' '.join(Vnames))
print("Initial date range:",mindate,"-",maxdate)
zconf=norm.ppf((1+conf)/2)
maxmult=20
numv=len(Vnames)
cogdate=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile+'.gz')).strftime('%Y-%m-%d')

def treeclassify(mutations):
  if '|synSNP:C14599T|' in mutations and '|synSNP:C3241T|' in mutations: return "XE"
  if '|S:R346T|' in mutations:
    if '|S:N658S|' in mutations: return "BA.4.6"
    elif '|S:Y248N|' in mutations: return "BA.2.76"
    elif '|S:L452M|' in mutations: return "BA.2.74"
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
          else:
            if '|S:A1020S|' in mutations and '|ORF7a:H47Y|' in mutations: return "BF.5"
            return "BA.5.2.1"
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
  if data[date][0]>=args.mincount and max(data[date][1:])>=args.mincount:
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

# Do simple regression to get decent initial values for Dirichlet-Multinomial regression

# Expand multinomial log probability about maximum:
# log(prod_r p_r ^ n_r) = constant - (1/2)*[ sum_r n_r.x_r^2 - (sum_r n_r.x_r)^2/(sum_r n_r) ] + O(x_^3)
# where x_r=log(p_r/n_r)
# Note that this is invariant under x_r -> x_r + constant, and we're expanding about the line (c,c,c,...,c).

# Let t = timestep, r = variant number
# Calculate, for each t, the quadratic form QF(x) = sum_r n_{t,r}x_{t,r}^2 - 1/(sum_r n_{t,r})*(sum_r n_{t,r}x_{t,r})^2
# Note that log likelihood(x) = (-1/2) * QF(x)
# I.e., find coefficients Q_{t,r,s} such that QF(x) = sum_{r,s} Q_{t,r,s}x_{t,r.x_{t,s}
Q=np.zeros([n,numv,numv])
Q[:,np.arange(numv),np.arange(numv)]=NN
Q-=NN[:,None,:]*NN[:,:,None]/NN.sum(axis=1)[:,None,None]

# We're looking for p_{t,r} = exp(a_r+t*b_r)/(constant independent of r), where b_r is the growth of variant r, so 
# substitute x_{t,r} = a_r + t*b_r - log(n_{t,r}), and sum over t.
# Amalgamate indexes corresponding to a_r and b_s as one single k index of size 2*numv.
# So (c_k) = (a_0,...,a_{numv-1},b_0,...,b_{numv-1}), and the quadratic (+linear) form, representing -2*LL, is
#    sum_{k,l} M_{k,l}.c_k.c_l - 2*sum_k R_k.c_k + constant, i.e. c^t.M.c - 2*R.c, with MLE(c) = M^{-1}R
M=np.zeros([numv*2,numv*2])
M[:numv,:numv]=Q.sum(axis=0)
M[numv:,:numv]=M[:numv,numv:]=(Q*np.arange(n)[:,None,None]).sum(axis=0)
M[numv:,numv:]=(Q*(np.arange(n)**2)[:,None,None]).sum(axis=0)

# Add gauge-fixing terms to M, requiring a_0=b_0=0 (other gauges are available).
scale=sqrt(NN.sum())# Some kind of size scale - result doesn't depend on this unless we change it by many orders of magnitude.
M[0,0]+=scale
M[numv,numv]+=scale
#M[:numv,:numv]+=scale*np.ones([numv,numv])
#M[numv:,numv:]+=scale*np.ones([numv,numv])

# Calculate linear term, R
LN=np.log(NN)
R=np.zeros(numv*2)
R[:numv]=(Q*LN[:,:,None]).sum(axis=(0,1))
R[numv:]=(Q*np.arange(n)[:,None,None]*LN[:,:,None]).sum(axis=(0,1))

# Calculate provisional covariance matrix and MLE
# Only differences of c[] elements and certain col,row differences of C[] make sense.
# E.g,. c[5]-c[2] and C[2,2]-C[2,5]-C[5,2]+C[5,5] make sense (assuming numv>5 in this example).
# Anything that can't be written as a function of these differences will be gauge-dependent.
C=np.linalg.inv(M)
c=C@R

# Calculate residuals (res), and overdispersion factor (mult), given n*(numv-1) degrees of freedom
# Effective covariance matrix, corrected for overdispersion, is then mult*C
res=c[None,:numv]+np.arange(n)[:,None]*c[None,numv:]-LN
mult=(Q*res[:,:,None]*res[:,None,:]).sum()/(n*(numv-1))
print("Residual(overdispersion) multiplier = %.3f"%mult,end="")
if mult<1: print("    - adjusting to 1");mult=1
else: print()
C*=mult

# xx[:numv] = intercepts
# xx[numv:2*numv] = growths
# xx[2*numv] = overdispersion multiplier
# NN[timestep up to n][variant up to numv] = count
def LL(xx,pr=0):
  LL=0
  aa=xx[:numv]
  bb=xx[numv:2*numv]
  mult=xx[2*numv]
  for t in range(n):
    rho=np.exp(aa+t*bb)
    N=NN[t].sum()
    s=max((N-mult)/(mult-1),1e-3)
    al=s*rho/rho.sum()
    DLL=gammaln(N+1)+gammaln(s)-gammaln(N+s)
    DLL+=gammaln(NN[t]+al).sum()-gammaln(al).sum()-gammaln(NN[t]+1).sum()
    if pr: print(mindate+t,NN[t],"logrho",np.log(rho),"al",al,"aa",aa,"bb",bb,"mult",mult,"DLL",DLL)
    LL+=DLL
  return LL

condition=1e3
def NLL(xx): return -LL(xx)/condition

def Hessian(xx,eps):
  N=len(xx)
  H=np.zeros([N,N])
  for i in range(N-1):
    for j in range(i+1,N):
      v=0
      for (s1,s2) in [(-1,-1),(-1,1),(1,-1),(1,1)]:
        x=np.copy(xx)
        x[i]+=s1*eps[i]
        x[j]+=s2*eps[j]
        v+=s1*s2*LL(x)
      H[i,j]=H[j,i]=v/(4*eps[i]*eps[j])
  for i in range(N):
    x=np.copy(xx)
    v=0
    for s in [-1,0,1]:
      x=np.copy(xx)
      x[i]+=s*eps[i]
      v+=(s*s*3-2)*LL(x)
    H[i,i]=v/eps[i]**2
  return H

if args.DM:
  bounds=[(c[i]-c[0]-3,c[i]-c[0]+3) for i in range(numv)]+[(c[i]-c[numv]-0.2,c[i]-c[numv]+0.2) for i in range(numv,2*numv)]+[(1.01,maxmult)]
  bounds[0]=bounds[numv]=(0,0)
  res=minimize(NLL,list(c)+[min(2,maxmult)],bounds=bounds, method="SLSQP", options={'ftol':1e-20, 'maxiter':10000})
  if not res.success: raise RuntimeError(res.message)
  print("Log likelihood: %.3f"%(LL(res.x)))
  xx=res.x
  for i in range(len(xx)):
    if bounds[i][0]<bounds[i][1] and (xx[i]<bounds[i][0]+1e-3 or xx[i]>bounds[i][1]-1e-3):
      if i<2*numv: desc=Vnames[i%numv]+[" intercept"," growth"][i//numv]
      else: desc="multiplier"
      print("Warning:",desc,"hit bound")
  print("Residual(overdispersion) multiplier from DM = %.3f"%xx[2*numv])
  
  eps=[1e-3]*numv+[1e-4]*numv+[1e-3]
  H=Hessian(xx,eps)
  # Add gauge-fixing terms (to ensure a0=b0=0). These won't affect the result after C is processed (below) to only see differences.
  H[0,0]+=scale
  H[numv,numv]+=scale
  C=-np.linalg.inv(H)
  c=xx

# Only gauge-invariant parts of C make sense: covariances referring to differences in intercepts and growth rates
# Now make a choice that all variances/covariances refer to baseline variant
# Convert query Cov[a_i,a_j] into Cov[a_i-a_0,a_j-a_0]
#               Cov[a_i,b_j] into Cov[a_i-a_0,b_j-b_0]
#               Cov[a_i,mult] into Cov[a_i-a_0,mult]
#               etc
C[:numv,      :]-=C[0,:][None,:]
C[numv:2*numv,:]-=C[numv,:][None,:]
C[:,      :numv]-=C[:,0][:,None]
C[:,numv:2*numv]-=C[:,numv][:,None]

# Sampling is most convenient way of getting CrIs for crossover points
numsamp=100000
test=mvn.rvs(mean=c,cov=C,size=numsamp)

out=[None]
for i in range(1,numv):
  print()
  grad=c[numv+i]
  graderr=sqrt(C[numv+i,numv+i])*zconf
  yoff=c[i]
  t=sorted(-test[:,i]/test[:,numv+i])
  cross=t[int(numsamp/2)]
  crosslow,crosshigh=t[int(numsamp*(1-conf)/2)],t[int(numsamp*(1+conf)/2)]
  crosserr=(crosshigh-crosslow)/2
  growthstr=f"Relative growth in {Vnames[i]} vs {Vnames[0]} of {grad*100:.1f}% ({(grad-graderr)*100:.1f}% - {(grad+graderr)*100:.1f}%) per day"
  doubstr=f"Doubling of ratio {Vnames[i]}/{Vnames[0]} every {log(2)/grad:.1f} ({log(2)/(grad+graderr):.1f} - {log(2)/(grad-graderr):.1f}) days"
  print(growthstr)
  print(doubstr)
  cr1=int(floor(cross))
  crossstr=f"Est'd %s/%s crossover on %s.%d +/- %.1f days"%(Vnames[i],Vnames[0],mindate+cr1,int((cross-cr1)*10),crosserr)
  print(crossstr)
  out.append((grad,graderr,yoff,cross,crosserr,growthstr,doubstr,crossstr))

datafn=location+'_%s'%('_'.join(Vnames))
future1=30
future2=30
maxt0=maxdate-mindate+1
maxt=maxt0+max(future1,future2)
samp=test[:,:numv,None]+test[:,numv:2*numv,None]*np.arange(maxt)[None,None,:]
b=test[:,numv:2*numv,None]
e=np.exp(samp)
eb=(e*b).sum(axis=1)/e.sum(axis=1)
beta_mean=np.mean(eb,axis=0)
zeropoint=beta_mean[maxt0-1]
beta_mean-=zeropoint
beta_low=np.quantile(eb,(1-conf)/2,0)-zeropoint
beta_high=np.quantile(eb,(1+conf)/2,0)-zeropoint
visthr=1e-6
ymin=50;ymax=-50
with open(datafn,'w') as fp:
  for t in range(maxt):
    if t<maxt0:
      v=VV[t]
      print(mindate+t,' '.join("%6d"%x for x in v),end='',file=fp)
    else:
      print(mindate+t,' '.join("%6s"%'-' for x in v),end='',file=fp)
    print(" %12g %12g %12g"%(beta_mean[t],beta_low[t],beta_high[t]),end='',file=fp)
    for i in range(numv):
      mu=c[i]+t*c[numv+i]
      var=C[i,i]+2*t*C[i,numv+i]+t**2*C[numv+i,numv+i]
      q0=mu-zconf*sqrt(var)
      q1=mu+zconf*sqrt(var)
      if t<maxt0:
        a,b=v[0]+1e-30,v[i]+1e-30
        y=log(b/a)
        prec=sqrt(a*b/(a+b))# Heuristic to guide size of blobs in plot
        print(" %12g %12g"%(y,prec),end='',file=fp)
        if prec>visthr:
          if args.plotpoints:
            ymin=min(ymin,y);ymax=max(ymax,y)
          if args.plotbands:
            ymin=min(ymin,q0)
            ymax=max(ymax,q1)
      else:
        print(" %12s %12s"%('-','-'),end='',file=fp)
      print(" %12g %12g"%(q0,q1),end='',file=fp)# Would be better as a formula in gnuplot, then it's continuous and can be extrapolated (though the code would be more cluttered)
      
    print(file=fp)
print("Written data to",datafn)

graphfn=datafn+'.png'
allothers=', '.join(Vnames[1:])
oneother=' or '.join(Vnames[1:])
if numv==2:
  number="it is"
  possessive="its"
else:
  number="they are"
  possessive="their"

graphtitle=f"New cases per day in the UK of {allothers} compared with {Vnames[0]}"
if future1>0: graphtitle+=f", with a {future1}-day projection"
graphtitle+=f"\\nNB: This is the est'd relative growth of {allothers} compared to {Vnames[0]}, not {possessive} absolute growth. It indicates how fast {number} taking over from {Vnames[0]}\\n"
if args.plotpoints:
  graphtitle+="Larger blobs indicate more certainty (more samples). "
graphtitle+=f"Description/caveats/current graph: http://sonorouschocolate.com/covid19/index.php/UK\\\\_variant\\\\_comparison\\nSource: Sequenced cases from COG-UK {cogdate}"
  
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
set title "{graphtitle}"
min(a,b)=(a<b)?a:b
"""
if future1>0 and not args.plotpoints:
  cmd+=f"""set arrow from "{maxdate}",graph 0 to "{maxdate}",graph 1 nohead lc 8 dashtype (40,20)\n"""
cmd+=f"""
plot [:"{str(maxdate+future1)}"] [{(ymin-0.1)/log(2)}:{max(ymax+(ymax-ymin)*(numv*0.1+0.1),0.5)/log(2)}]"""
#plot [:] [{(ymin-0.5)/log(2)}:{max(ymax+0.8*numv-0.6,1.8)/log(2)}]"""

for i in range(1,numv):
  (grad,graderr,yoff,cross,crosserr,growthstr,doubstr,crossstr)=out[i]
  if args.plotpoints: cmd+=f""" "{datafn}" u 1:((${numv+5+4*i})/log(2)):(min(${numv+5+4*i+1},20)/{maxt0/8.}) pt 5 lc {i} ps variable title "","""
  if args.plotbands: cmd+=f""" "{datafn}" u 1:((${numv+5+4*i+2})/log(2)):((${numv+5+4*i+3})/log(2)) w filledcurves lc {i} title "","""
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
  graderr=sqrt(C[numv+i,numv+i])*zconf
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


graphfn=datafn+".growthproj.png"
linetitle=f"The part of the change (compared with {maxdate}) in the percentage daily increase\\n in new cases per day that is due to the change in variant mixture"
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
set ylabel "Percentage growth rate of new cases per day, compared with {maxdate}"
set style fill transparent solid 0.25
set style fill noborder
set format y "%.1f%%"

set output "{graphfn}"
set title "Estimated effect of variant mixture {', '.join(Vnames)} on the overall growth rate in new cases/day, set at 0 on {maxdate}\\nNB: This growth rate is affected by other things - only the contribution to the growth rate due to the variant mixture is shown here\\n"""
cmd+=f"""Description/caveats/current graph: http://sonorouschocolate.com/covid19/index.php/UK\\\\_variant\\\\_comparison\\nSource: Sequenced cases from COG-UK {cogdate}"
set arrow from "{maxdate}",graph 0 to "{maxdate}",graph 1 nohead lc 8 dashtype (40,20)
plot [:"{str(maxdate+future2)}"] """
cmd+=f""" "{datafn}" u 1:((${numv+2})*100) lc 1 lw 2 w lines title "{linetitle}", """
cmd+=f""" "{datafn}" u 1:((${numv+3})*100):((${numv+4})*100) lc 1 w filledcurves title "" """

po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE)
p=po.stdin
p.write(cmd.encode('utf-8'))
p.close()
po.wait()
print()
print("Written graph to",graphfn)
print()

