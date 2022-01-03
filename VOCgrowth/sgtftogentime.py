from stuff import *
from math import log,exp
from scipy.optimize import minimize
from random import randrange
from math import sqrt
import numpy as np
np.set_printoptions(precision=6,linewidth=250,suppress=True)

minday=Date('2021-11-25')
maxday=Date('2021-12-25')# Only go up to dates strictly before this one
pubday=getpublishdate()
discard=3# Discard last few case counts by specimen date since these are incomplete (irrelevant here because we're stopping much earlier anyway)
mincount=10
step=7
outdir='gentimeoutput'
os.makedirs(outdir,exist_ok=True)

l=[x for x in os.listdir('.') if x[:19]=='sgtf_regionepicurve']
if l==[]: raise RuntimeError("No sgtf_regionepicurve*.csv file found in current directory; download from https://www.gov.uk/government/publications/covid-19-omicron-daily-overview")
sgtf=loadcsv(max(l))
regions=sorted(list(set(sgtf['UKHSA_region'])))

ages=[(0,150)]
casesbyregion={}
for loc in regions:
  sp0,sp=getcasesbyagespeccomplete(minday=minday,maxday=pubday,ages=ages,location=loc)
  casesbyregion[loc]=sp[:pubday-discard-minday]
casesbyregion['England']=sum(casesbyregion.values())
nspec=casesbyregion['England'].shape[0]
  
# Get SGTF data into a suitable form
vocnum={}
background=[0,0]
nsgtf=max(Date(d) for d in sgtf['specimen_date'])+1-minday
for (date,region,variant,n) in zip(sgtf['specimen_date'],sgtf['UKHSA_region'],sgtf['sgtf'],sgtf['n']):
  day=Date(date)
  daynum=day-minday
  if daynum>=0 and daynum<nsgtf:
    place=region
    if place not in vocnum: vocnum[place]=np.zeros([nsgtf,2],dtype=int)
    vocnum[place][daynum][int("SGTF" in variant)]+=n
  if date>=Date('2021-10-01') and date<Date('2021-11-10'):
    background[int("SGTF" in variant)]+=n

# Adjust for non-Omicron SGTFs, based on the assumption that these are in a non-location-dependent proportion to the number of non-Omicron cases
f=background[1]/background[0]
for place in vocnum:
  for daynum in range(nsgtf):
    vocnum[place][daynum][1]=max(vocnum[place][daynum][1]-int(f*vocnum[place][daynum][0]+.5),0)
vocnum['England']=sum(vocnum.values())

def regressYonX(W,X,Y):
  W=np.array(W);X=np.array(X);Y=np.array(Y)
  m=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
  r=np.array([sum(W*Y),sum(W*X*Y)])
  return np.linalg.solve(m,r)
  

date0=Date('2021-11-25')
date1=Date('2021-12-29')
vn=vocnum['England']
dat=[]
for date in Daterange(date0,date1):
  vv=vn[date-minday]
  if vv[0]>0 and vv[1]>0:
    dat.append((str(date),log(vv[1]/vv[0]),1/(1/vv[0]+1/vv[1])))
fn=os.path.join(outdir,'deltaomicron.dat')
with open(fn,'w') as fp:
  for (date,y,w) in dat:
    print(date,'%10.6f %12g'%(y,w),file=fp)

cutoff=Date('2021-12-14')
ll=[]
for after in [0,1]:
  W=[];X=[];Y=[]
  for (x,(date,y,w)) in enumerate(dat):
    if (date>cutoff)==after: W.append(w);X.append(x);Y.append(y)
  ll.append(regressYonX(W,X,Y))

cmd="""
set terminal pngcairo font "sans,13" size 1920,1920
set bmargin 5;set lmargin 15;set rmargin 10;set tmargin 5
cd "%s"
set output "deltaomicron.png"
set key left Left reverse
set xdata time;fmt="%%Y-%%m-%%d";set timefmt fmt;set format x "%%Y-%%m-%%d"
set title "Simple regression of log(Omicron)/log(Delta) against time\\nDiscussion: http://sonorouschocolate.com/covid19/index.php?title=Estimating\\\\_Generation\\\\_Time\\\\_Of\\\\_Omicron\\nData source: SGTF numbers from UKHSA Omicron daily overview"
set ylabel "log(number of likely Omicron on given day/number of likely Delta on given day)"
set xtics "2020-01-06", 604800
plot "deltaomicron.dat" u 1:2:(sqrt($3)/15) w points pt 5 ps variable lc 2 title "log(number of likely Omicron on given day/number of likely Delta on given day); larger sizes indicates more confidence in vertical position",%g+%g*(x-strptime(fmt,"%s"))/86400 lc 3 lw 3 title "Best fit line up to %s: slope %.3f",%g+%g*(x-strptime(fmt,"%s"))/86400 lc 4 lw 3 title "Best fit line after %s: slope %.3f"
"""%(outdir,ll[0][0],ll[0][1],str(minday),str(cutoff),ll[0][1],ll[1][0],ll[1][1],str(minday),str(cutoff),ll[1][1])
po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE)
p=po.stdin
p.write(cmd.encode('utf-8'))
p.close()
po.wait()

data=[]# list of (growth in Delta, growth in Omicron)
var=[]# list of variances of the above
with open('gg-by-region%d'%step,'w') as fp:
  n=maxday-minday
  #n=min(nsgtf,nspec)
  for loc in vocnum:
    if loc=='England': continue
    vv=vocnum[loc][:n,:]
    cv=casesbyregion[loc][:n,:].sum(axis=1)[:,None]*(vv/vv.sum(axis=1)[:,None])
    #cv=vv
    for i in range(n-step):
      if (vv[i:i+2*step:step,:]>=mincount).all():
        gr0=log(cv[i+step][0]/cv[i][0])/step
        gr1=log(cv[i+step][1]/cv[i][1])/step
        v0=1/vv[i,0]+1/vv[i+step,0]
        v1=1/vv[i,1]+1/vv[i+step,1]
        data.append((gr0,gr1))
        var.append((v0,v1))
        print("%7.4f %7.4f"%(data[-1]),Date(minday+i),Date(minday+i+step),"%6d %6d %6d %6d"%tuple(vv[i:i+2*step:step,:].reshape(-1)),loc,file=fp)
data=np.array(data)
var=np.array(var)

# 2d weighted regression
def regress(data,var):
  # To get a starting point, do standard weighted regression, pretending uncertainty in (x,y) is all in the y
  W=1/(var.sum(axis=1))
  (X,Y)=data[:,0],data[:,1]
  c=regressYonX(W,X,Y)

  # Calculate weighted square-distance to putative best fit line of gradient b
  def f(b):
    (X,Y)=data[:,0],data[:,1]
    (V,W)=var[:,0],var[:,1]
    N=Y-b*X
    D=W+b*b*V
    a=sum(N/D)/sum(1/D)
    return (a,sum((N-a)*(N-a)/(2*D)))
  
  f0=f(c[1])[1]
  def g(b): return f(b)[1]/f0
  
  bounds=(c[1]*0.5,c[1]/0.5)
  res=minimize(g,[c[1]],bounds=[bounds],method="SLSQP")
  if not res.success: raise RuntimeError(res.message)
  b=res.x[0]
  (a,s)=f(b)
  #if b<bounds[0]+1e-6 or b>bounds[1]-1e-6: raise RuntimeError("Hit bounds %g %g"%bounds)
  #if b<bounds[0]+1e-6 or b>bounds[1]-1e-6: print("Warning: Hit bounds %g %g"%bounds,file=sys.stderr)
  return (a,b)

# Bootstrap to get confidence intervals
n=len(data)
nsamp=1000
conf=0.95

samples=[]
print("Generating %d samples"%nsamp)
for samp in range(nsamp):
  data1=np.zeros(data.shape)
  var1=np.zeros(var.shape)
  for i in range(n):
    r=randrange(n)
    data1[i]=data[r]
    var1[i]=var[r]
  samples.append(regress(data1,var1))

# Compare with Delta

# https://www.medrxiv.org/content/10.1101/2021.10.21.21265216v1
# Delta (intrinsic): 4.6 (4.0-5.4) days  3.1 (3.0-3.7) days
# >>> from scipy.stats import gamma
# >>> a=170;mu=4.64;gamma.ppf(.025,a,scale=mu/a),gamma.ppf(.975,a,scale=mu/a)
# (3.9687044779630987, 5.3629770162251926)

from scipy.stats import gamma,expon
csamples=[]
for (a,b) in samples:
  # Generate meanGT_Delta with characteristics 4.6 (4.0 - 5.4)
  ag=170;mu=4.64;mean_d=gamma.rvs(ag,scale=mu/ag)

  # Generate sdGT_Delta with characteristics 3.1 (3.0 - 3.7)
  # This is very skew, so do it in two pieces, setting the median (not the mean) to 3.1
  e=expon.ppf(0.95)
  if randrange(2): sd_d=3.1+.6/e*expon.rvs()
  else: sd_d=3.1-.1/e*expon.rvs()
  
  th_d=sd_d**2/mean_d
  k_d=(mean_d/sd_d)**2
  th_o=1/(b/th_d-a)
  k_o=k_d# Assume equal shape (or coefficient of variation)
  rr=(b*th_o/th_d)**k_d
  mean_o=th_o*k_o
  sd_o=th_o*sqrt(k_o)
  gtr=mean_o/mean_d
  csamples.append([mean_d,sd_d,mean_o,sd_o,gtr,rr])
csamples=np.array(csamples)

desc=["meanGT_Delta","sdGT_Delta","meanGT_Omicron","sdGT_Omicron","GT_Om/GT_Delta","R_Om/R_Delta"]
nd=len(desc)
for i in range(nd):
  lsamp=list(csamples[:,i])
  lsamp.sort()
  low=lsamp[int((1-conf)/2*nsamp)]
  med=lsamp[int(0.5*nsamp)]
  high=lsamp[int((1+conf)/2*nsamp)]
  print("%-16s  %6.3f (%.3f - %.3f)"%(desc[i],med,low,high))
  
# Simple CI of gradient
samples.sort(key=lambda x:x[1])
low=samples[int((1-conf)/2*nsamp)]
high=samples[int((1+conf)/2*nsamp)]
central=regress(data,var)
print("Central estimate: y=%.3f+%.3f*x (y=%.3f+%.3f*x - y=%.3f+%.3f*x)"%(central+low+high))
(a,b)=central

# Make lower, upper curves for plotting
x0,x1=data[:,0].min(),data[:,0].max()
d=x1-x0
x0-=0.05*d
x1+=0.05*d
n=1000
probsamples=np.array(samples[int((1-conf)/2*nsamp):int((1+conf)/2*nsamp)])
low=[];high=[]
xrange=np.arange(x0,x1,(x1-x0)/n)
CI=[]
for x in xrange:
  yy=probsamples@[1,x]
  CI.append([x,yy.min(),yy.max()])

with open(os.path.join(outdir,'rawdata'),'w') as fp:
  W=1/(var.sum(axis=1))
  for ((x,y),w) in zip(data,W): print("%10.7f %10.7f    %10.3f"%(x,y,w),file=fp)
  
with open(os.path.join(outdir,'CI'),'w') as fp:
  for (x,y0,y1) in CI:
    print("%10.7f %10.7f %10.7f"%(x,y0,y1),file=fp)

(a,b)=central
cmd="""
set terminal pngcairo font "sans,13" size 1920,1920
set bmargin 5;set lmargin 15;set rmargin 10;set tmargin 5
cd "%s"
set output "growthcomparison.png"
set key left
set title "Comparison of growth rate of Omicron and growth rate of Delta\\nDiscussion: http://sonorouschocolate.com/covid19/index.php?title=Estimating\\\\_Generation\\\\_Time\\\\_Of\\\\_Omicron\\nData sources: UKHSA Omicron daily overview and dashboard"
set xlabel "Average growth per day of Delta over 7-day period"
set ylabel "Average growth per day of Omicron over 7-day period"
set style fill transparent solid 0.5 noborder
plot "rawdata" u 1:2:(sqrt($3)/10) w points pt 5 ps variable lc 2 title "Pairs of growths over 7 day intervals over regions of England and days in December; larger sizes indicate more confidence in position", "CI" u 1:2:3 w filledcurves lc 3 title "95%% bootstrap confidence interval of best fit line",%g+%g*x lc 3 lw 3 title "Best fit line: y=%.3f+%.3f*x"
"""%(outdir,a,b,a,b)
po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE)
p=po.stdin
p.write(cmd.encode('utf-8'))
p.close()
po.wait()
