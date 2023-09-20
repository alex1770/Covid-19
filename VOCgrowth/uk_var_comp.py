import sys,pickle,os,argparse
from stuff import *
from scipy.stats import norm, multivariate_normal as mvn
from scipy.optimize import minimize
from scipy.special import gammaln,digamma
import numpy as np
from math import sqrt,floor,log,exp
from variantaliases import aliases
from classify import classify, contractlin, expandlin
import hashlib

np.set_printoptions(precision=6,suppress=True,linewidth=200)

cachedir="cogukcachedir"
datafile='cog_metadata.csv'
conf=0.95

parser=argparse.ArgumentParser()
parser.add_argument('-c',  '--mincount',    type=int,default=0,    help="Minimum variant count considered")
parser.add_argument('-d',  '--decluster',   action="store_true",   help="Decluster the variant counts to reduce the effect of targetted sequencing")
parser.add_argument('-s',  '--simple',      action="store_true",   help="Use simple weighted and overdispersion-corrected logistic regression instead of Dirichlet-Multinomial")
parser.add_argument('-f',  '--mindate',     default="2022-01-01",  help="Min sample date of sequence")
#parser.add_argument('-g',  '--gisaid',      action="store_true",   help="Use GISAID input")
parser.add_argument('-t',  '--maxdate',     default="9999-12-31",  help="Max sample date of sequence")
parser.add_argument('-l',  '--lineages',    default="BA.4*,BA.5*", help="Comma-separated list of lineages/variants")
parser.add_argument('-m',  '--maxmult',     type=float,default=20, help="Maximum overdispersion multiplier")
parser.add_argument('-p',  '--plotpoints',  action="store_true",   help="Whether to plot points corresponding to log(num(variant)/num(base variant))")
parser.add_argument('-b',  '--plotbands',   action="store_true",   help="Whether to plot confidence bands around best-fit lines")
parser.add_argument('-F',  '--future',      type=int,default=30,   help="Number of days ahead to project")
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
numv=len(Vnames)
cogdate=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile+'.gz')).strftime('%Y-%m-%d')

# lin is assumed to have already been expanded
def patmatch(lin):
  ind=-1
  for i in range(len(Vnames)):
    exact=targlinsexact[i]
    prefix=targlinsprefix[i]
    if lin==exact:
      ind=i
      if prefix=='-': return ind# Exact match with non-wildcard takes precedence over anything later
    if lin[:len(prefix)]==prefix: ind=i
  return ind

jn='_'.join(Vnames)
if len(jn)>200: jn=hashlib.sha256(jn.encode("utf-8")).hexdigest()[-16:]
if args.decluster: jn+="__decluster"
fn=os.path.join(cachedir,location+'_'+jn+'__%s'%cogdate)
if os.path.isfile(fn):
  with open(fn,'rb') as fp:
    counts=pickle.load(fp)
else:
  # Wildcard ending is replaced with '.'. It's a match if it's equal to a prefix of (database lineage)+'.'
  # Note that BA.5* will match BA.5 and BA.5.1 but not BA.53
  #           BA.5.* will match BA.5.1 but not BA.5 or BA.53
  targlinsexact=[]
  targlinsprefix=[]
  for lin0 in Vnames:
    lin=expandlin(lin0)
    if lin=='*' or lin[-2:]=='.*': exact="-";prefix=lin[:-1]
    elif lin[-1]=='*': exact=lin[:-1];prefix=lin[:-1]+'.'
    else: exact=lin;prefix="-"
    targlinsexact.append(expandlin(exact))
    targlinsprefix.append(expandlin(prefix))
  counts={}
  seqdate=None
  for (name,date,p2,lin,mutations,ambiguities) in csvrows(datafile,['sequence_name','sample_date','is_pillar_2','usher_lineage','mutations','ambiguities']):
    #if p2!='Y': continue
    if mutations=="": continue
    if location!="UK":
      country=name.split('/')[0]
      if country!=location: continue
    if not (len(date)==10 and date[:2]=="20" and date[4]=="-" and date[7]=="-"): continue
    
    # Try to assign sublineage to one of the given lineages. E.g., if Vnames=["BA.1*","BA.1.1*","BA.2"] then BA.1.14 is counted as BA.1* but BA.1.1.14 is counted as BA.1.1*
    ind=patmatch(expandlin(lin))
    if ind==-1: continue
    if args.decluster:
      if date!=seqdate:
        if seqdate!=None: counts[seqdate]=[declusternumber(s) for s in seqlist]
        seqlist=[[] for _ in range(numv)]
        seqdate=date
        print(date)
      seqlist[ind].append(conv_climb_metadata_to_sequence(mutations,ambiguities))
    else:
      if date not in counts: counts[date]=[0]*numv
      counts[date][ind]+=1
    #if date<mindate: break# If we don't abort early, then can store a cache file that works for any date range - alter
  if args.decluster:
    if seqdate!=None: counts[seqdate]=[declusternumber(s) for s in seqlist]
  os.makedirs(cachedir,exist_ok=True)
  with open(fn,'wb') as fp:
    pickle.dump(counts,fp)

mindate1=Date('2099-12-31')
maxdate1=Date('2000-01-01')
for date in counts:
  if counts[date][0]>=args.mincount and max(counts[date][1:])>=args.mincount:
    mindate1=min(mindate1,Date(date))
    maxdate1=max(maxdate1,Date(date))
mindate=max(mindate,mindate1)
maxdate=min(maxdate,maxdate1)
print("Reduced date range:",mindate,"-",maxdate)

VV=[]
for date in Daterange(mindate,maxdate+1):
  v=counts.get(str(date),[0]*numv)
  VV.append(v)
if VV==[]: print("No data points found");sys.exit(0)
VV=np.array(VV)
npd=(VV>0).sum(axis=0);bad=0
for i in range(numv):
  if npd[i]<2: print("Variant",Vnames[i],"only has positive counts on",npd[i],"day"+"s"*int(npd[i]!=1));bad=1
if bad: raise RuntimeError("Variant counts too low")

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
# I.e., find coefficients Q_{t,r,s} such that QF(x) = sum_{r,s} Q_{t,r,s}x_{t,r}x_{t,s}
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

def dLL(xx,pr=0):
  dLL=0
  aa=xx[:numv]
  bb=xx[numv:2*numv]
  mult=xx[2*numv]
  for t in range(n):
    rho=np.exp(aa+t*bb)
    N=NN[t].sum()
    s=max((N-mult)/(mult-1),1e-3)
    al=s*rho/rho.sum()
    dal_da=s*(rho.sum()*np.diag(rho)-np.outer(rho,rho))/rho.sum()**2
    dal_db=dal_da*t
    if s<1e-3: dal_dm=np.zeros(numv)[:,None]
    else: dal_dm=(rho/rho.sum()*-(N-1)/(mult-1)**2)[:,None]
    dal=np.concatenate([dal_da,dal_db,dal_dm],axis=1)
    ds=np.array([0]*(numv*2)+[-(N-1)/(mult-1)**2])
    if pr: print(mindate+t,((digamma(s)-digamma(N+s))*ds+(digamma(NN[t]+al)-digamma(al))@dal)[2])
    dLL+=(digamma(s)-digamma(N+s))*ds+(digamma(NN[t]+al)-digamma(al))@dal
  return dLL

condition=1e3
def NLL(xx): return -LL(xx)/condition
def NdLL(xx): return -dLL(xx)/condition

def Hessian(xx,eps):
  N=len(xx)
  H=np.zeros([N,N])
  for i in range(N):
    leps=[0]*i+[eps[i]]+[0]*(N-1-i)
    H[i,:]=(dLL(xx+leps)-dLL(xx-leps))/(2*eps[i])
  return (H+H.T)/2

if not args.simple:
  bounds=[(c[i]-c[0]-30,c[i]-c[0]+30) for i in range(numv)]+[(c[i]-c[numv]-0.2,c[i]-c[numv]+0.2) for i in range(numv,2*numv)]+[(min(1.01,args.maxmult),args.maxmult)]
  bounds[0]=bounds[numv]=(0,0)
  res=minimize(NLL,list(c)+[min(2,args.maxmult)],bounds=bounds, jac=NdLL, method="SLSQP", options={'ftol':1e-10, 'maxiter':10000})
  if not res.success: raise RuntimeError(res.message)
  print("Log likelihood: %.3f"%(LL(res.x)))
  xx=res.x
  odmbound=False
  for i in range(len(xx)):
    if bounds[i][0]<bounds[i][1] and (xx[i]<bounds[i][0]+1e-3 or xx[i]>bounds[i][1]-1e-3):
      if i<2*numv: print("Error:",Vnames[i%numv]+[" intercept"," growth"][i//numv],"hit bound")
      else: odmbound=True;print("Note: overdispersion multiplier hit bound")
  print("Residual(overdispersion) multiplier from DM = %.3f"%xx[2*numv])
  
  eps=np.array([1e-5]*numv+[1e-6]*numv+[1e-5])
  H=Hessian(xx,eps)
  # Add gauge-fixing terms (to ensure a0=b0=0). These won't affect the result after C is processed (below) to only see differences.
  H[0,0]-=scale
  H[numv,numv]-=scale
  if odmbound:# If multiplier hit bound, then need to treat it as fixed
    C=np.zeros([2*numv+1,2*numv+1])
    C[:2*numv,:2*numv]=-np.linalg.inv(H[:2*numv,:2*numv])
    C[2*numv,2*numv]=0
  else:
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
numsamp=50000
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

datafn=location+'_%s'%jn
maxt0=maxdate-mindate+1
maxt=maxt0+args.future
samp=test[:,:numv,None]+test[:,numv:2*numv,None]*np.arange(maxt)[None,None,:]
b=test[:,numv:2*numv,None]
e=np.exp(samp)
eb=(e*b).sum(axis=1)/e.sum(axis=1)
eb2=(e*b*b).sum(axis=1)/e.sum(axis=1)
vb=eb2-eb*eb
zeropoint=eb[:,maxt0-1]
growth_mean=np.mean(eb-zeropoint[:,None],axis=0)
growth_low=np.quantile(eb-zeropoint[:,None],(1-conf)/2,0)
growth_high=np.quantile(eb-zeropoint[:,None],(1+conf)/2,0)
pressure_mean=np.mean(vb,axis=0)
pressure_low=np.quantile(vb,(1-conf)/2,0)
pressure_high=np.quantile(vb,(1+conf)/2,0)
visthr=1e-6
ymin=50;ymax=-50
with open(datafn,'w') as fp:
  for t in range(maxt):
    if t<maxt0:
      v=VV[t]
      print(mindate+t,' '.join("%8g"%x for x in v),end='',file=fp)
    else:
      print(mindate+t,' '.join("%8s"%'-' for x in v),end='',file=fp)
    print(" %12g %12g %12g"%(growth_mean[t],growth_low[t],growth_high[t]),end='',file=fp)
    print(" %12g %12g %12g"%(pressure_mean[t],pressure_low[t],pressure_high[t]),end='',file=fp)
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
if args.future>0: graphtitle+=f", with a {args.future}-day projection"
graphtitle+=f"\\nNB: This is the est'd relative growth of {allothers} compared to {Vnames[0]}, not {possessive} absolute growth. It indicates how fast {number} taking over from {Vnames[0]}\\n"
if args.plotpoints:
  graphtitle+="Larger blobs indicate more certainty (more samples). "
graphtitle+=f"Description/caveats/current graph: http://sonorouschocolate.com/covid19/index.php/UK\\\\_variant\\\\_comparison\\nSource: Sequenced cases from CLIMB {cogdate}"
  
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
if args.future>0:# and not args.plotpoints:
  cmd+=f"""set arrow from "{maxdate}",graph 0 to "{maxdate}",graph 1 nohead lc 8 dashtype (40,20)\n"""
cmd+=f"""
plot [:"{str(maxdate+args.future)}"] [{(ymin-0.1)/log(2)}:{max(ymax+(ymax-ymin)*(numv*0.1+0.1),0.5)/log(2)}]"""
#plot [:] [{(ymin-0.5)/log(2)}:{max(ymax+0.8*numv-0.6,1.8)/log(2)}]"""

for i in range(1,numv):
  (grad,graderr,yoff,cross,crosserr,growthstr,doubstr,crossstr)=out[i]
  if args.plotpoints: cmd+=f""" "{datafn}" u 1:((${numv+8+4*i})/log(2)):(min(${numv+8+4*i+1},20)/{maxt0/8.}) pt 5 lc {i} ps variable title "","""
  if args.plotbands: cmd+=f""" "{datafn}" u 1:((${numv+8+4*i+2})/log(2)):((${numv+8+4*i+3})/log(2)) w filledcurves lc {i} title "","""
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

#for date in Daterange(mindate,maxdate+50):
#  print(date,end="")
#  q=[]
#  for i in range(numv):
#    if i>0: (grad,graderr,yoff,cross,crosserr,growthstr,doubstr,crossstr)=out[i]
#    else: grad=yoff=0
#    q.append(exp((date-mindate)*grad+yoff))
#  for i in range(numv):
#    print(" %6.1f%%"%(q[i]/sum(q)*100),end="")
#  print()
#print()

last=14
proj=42
AR=NN.sum(axis=0)/NN.sum()
PR=NN[-last:,:].sum(axis=0)/NN[-last:,:].sum()
stats=[]
for i in range(numv):
  grad=c[numv+i]
  graderr=sqrt(C[numv+i,numv+i])*zconf
  stats.append([i,grad,graderr,AR[i],PR[i],PR[i]*exp(grad*proj)])
stats=np.array(stats)
stats[:,5]=stats[:,5]/stats[:,5].sum()
orders=stats.argsort(axis=0)
types=[("original given order",0),("growth rate relative to %s"%Vnames[0],1),("relative prevalence over last %d days"%last,4),("relative prevalence projected forward %d days"%proj,5)]
for desc,col in types:
  print("Ordered by",desc)
  print("        Lineage ---------Growth---------     Alldays  Last%02ddays Proj%02ddays"%(last,proj))
  for i in range(numv):
    row=stats[orders[i,col],:]
    print("%15s %6.3f (%6.3f - %6.3f)      %5.1f%%      %5.1f%%     %5.1f%%"%(Vnames[int(row[0])],row[1],row[1]-row[2],row[1]+row[2],row[3]*100,row[4]*100,row[5]*100))
  print()


graphfn=datafn+".growthproj.png"
linetitle=f"Variant effect (integral of variant pressure; effect of changing variant mixture on the growth rate in new cases per day, compared with {maxdate})"
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
set title "Estimated effect of the changing variant mixture on the overall growth rate in new cases per day\\nVariants considered: {{/=10{', '.join(Vnames)}}}\\nNB: changes in growth rate can arise from several causes - only the contribution to the change in growth rate due to the variant mixture is shown here\\n"""
cmd+=f"""Description/caveats/current graph: http://sonorouschocolate.com/covid19/index.php/UK\\\\_variant\\\\_comparison\\nSource: Sequenced cases from CLIMB {cogdate}"
set arrow from "{maxdate}",graph 0 to "{maxdate}",graph 1 nohead lc 8 dashtype (40,20)
plot [:"{str(maxdate+args.future)}"] """
cmd+=f""" "{datafn}" u 1:((${numv+2})*100) lc 1 lw 2 w lines title "{linetitle}", """
cmd+=f""" "{datafn}" u 1:((${numv+3})*100):((${numv+4})*100) lc 1 w filledcurves title "" """

po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE)
p=po.stdin
p.write(cmd.encode('utf-8'))
p.close()
po.wait()
print()
print("Written graph to",graphfn)


graphfn=datafn+".variantpressure.png"
linetitle=f"Variant pressure (effect of changing variant mixture on the rate of change of the growth in overall new cases per day)"
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
set ylabel "Variant pressure (percentage growth rate per day per day)"
set style fill transparent solid 0.25
set style fill noborder
set format y "%.2f%%"

set output "{graphfn}"
set title "Estimated pressure the changing variant mixture puts on the change in overall growth rate of new cases per day, in new cases per day per day\\nVariants considered: {{/=10{', '.join(Vnames)}}}\\nNB: changes in growth rate can arise from several causes - only the contribution to the change in growth rate due to the variant mixture is shown here\\n"""
cmd+=f"""Description/caveats/current graph: http://sonorouschocolate.com/covid19/index.php/UK\\\\_variant\\\\_comparison\\nSource: Sequenced cases from CLIMB {cogdate}"
set arrow from "{maxdate}",graph 0 to "{maxdate}",graph 1 nohead lc 8 dashtype (40,20)
plot [:"{str(maxdate+args.future)}"] [:{min(max(pressure_mean)*2,max(pressure_high)*1.05)*100}] """
cmd+=f""" "{datafn}" u 1:((${numv+5})*100) lc 1 lw 2 w lines title "{linetitle}", """
cmd+=f""" "{datafn}" u 1:((${numv+6})*100):((${numv+7})*100) lc 1 w filledcurves title "" """

po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE)
p=po.stdin
p.write(cmd.encode('utf-8'))
p.close()
po.wait()
print()
print("Written graph to",graphfn)
print()
