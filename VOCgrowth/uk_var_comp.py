import sys
from stuff import *
from scipy.stats import norm

datafile='cog_metadata.csv'
VV=['BA.1','BA.2']
mindate=Date('2000-01-01')
maxdate=Date('2099-12-31')
mincount=5
conf=0.95

if len(sys.argv)>1: VV[0]=sys.argv[1]
if len(sys.argv)>2: VV[1]=sys.argv[2]
if len(sys.argv)>3: mindate=Date(sys.argv[3])
if len(sys.argv)>4: maxdate=Date(sys.argv[4])
zconf=norm.ppf((1+conf)/2)

data={}
for (date,p2,lin) in csvrows(datafile,['sample_date','is_pillar_2','lineage']):
  if lin not in VV or len(date)!=10: continue
  #if p2!='Y': continue
  i=VV.index(lin)
  if date not in data: data[date]=[0,0]
  data[date][i]+=1
  #if date<mindate: break

mindate1=Date('2099-12-31')
maxdate1=Date('2000-01-01')
for date in data:
  if data[date][0]>=mincount and data[date][1]>=mincount:
    mindate1=min(mindate1,Date(date))
    maxdate1=max(maxdate1,Date(date))
mindate=max(mindate,mindate1)
maxdate=min(maxdate,maxdate1)

datafn='UK_%s_%s'%tuple(VV)
V0=[];V1=[]
with open(datafn,'w') as fp:
  for date in Daterange(mindate,maxdate+1):
    v=data.get(str(date),[0,0])
    V0.append(v[0])
    V1.append(v[1])
    print(date,"%6d %6d"%tuple(v),file=fp)
print("Written data to",datafn)

import numpy as np
from math import sqrt,floor,log

# L(a,b) = -(1/2)sum_i w_i(a+bx_i-y_i)^2
#      c = (a,b); can write L(c)
# L is quadratic in a,b so
# dL/dc = dL/dc(0) + (d^2L/dc^2).c = r - M.c, where M and r are as below (known and constant, i.e. independent of a,b).
# So c*=M^{-1}r and L(x) = L(c*) - (1/2)(x-c*)^t.M.(x-c*)
# To correct for dependence, we're going to deem the average residual to be equal to 1, which we're going to achieve by rescaling V0, V1.
# So imagine: V0/=mult, V1/=mult, W/=mult, M/=mult, r/=mult, C*=mult, X, Y, a, b, c unchanged.
# Then M/mult is the precision (inverse-covariance) matrix = observed Fisher information,
# and (M/mult)^{-1} is the "observed covariance" matrix, and bottom right of this is est variance of b, the gradient.
# Effectively our posterior is (a,b) ~ MVN(c,C)
# We're interested in b (growth) and a/b (essentially the crossover point). Can get a/b by simulation
# or can get it simply and analytically if we don't worry about the small chance of b going negative.
V0=np.array(V0)+1e-30
V1=np.array(V1)+1e-30
n=len(V0)
W=V0*V1/(V0+V1)
X=np.arange(n)
Y=np.log(V1/V0)
M=np.array([[sum(W), sum(W*X)], [sum(W*X), sum(W*X*X)]])
r=np.array([sum(W*Y),sum(W*X*Y)])
C=np.linalg.inv(M)
c=C@r
res=c[0]+c[1]*X-Y
mult=(W*res*res).sum()/n
print("Residual multiplier = %.3f"%mult)
qa=c[1]**2-zconf**2*C[1,1]*mult
qb=c[0]*c[1]-zconf**2*C[0,1]*mult
qc=c[0]**2-zconf**2*C[0,0]*mult
descrim=qb**2-qa*qc
lam0=(-qb-sqrt(descrim))/qa
lam1=(-qb+sqrt(descrim))/qa
# from scipy.stats import multivariate_normal as mvn
# N=100000
# test=mvn.rvs(mean=c,cov=C,size=N)
# t=sorted(-test[:,0]/test[:,1])
# print(t[int(N*(1-conf)/2)],t[int(N*(1+conf)/2)])

grad=c[1]
graderr=sqrt(mult*C[1,1])*zconf
cross=-qb/qa
crosserr=sqrt(descrim)/qa
growthstr=f"Relative growth in {VV[1]} vs {VV[0]} of {grad:.3f} ({grad-graderr:.3f} - {grad+graderr:.3f}) per day"
doubstr=f"Doubling of {VV[1]}/{VV[0]} every {log(2)/grad:.1f} ({log(2)/(grad+graderr):.1f} - {log(2)/(grad-graderr):.1f}) days"
print(growthstr)
print(doubstr)
cr1=int(floor(cross))
crossstr="Crossover on %s.%d +/- %.1f days"%(mindate+cr1,int((cross-cr1)*10),crosserr)
print(crossstr)

cogdate=datetime.datetime.utcfromtimestamp(os.path.getmtime(datafile+'.gz')).strftime('%Y-%m-%d')
graphfn=datafn+'.png'

cmd=f"""
set xdata time
set key left Left reverse
set key spacing 3
fmt="%Y-%m-%d"
set timefmt fmt
set format x fmt
set xtics "2020-01-06", 604800
set xtics rotate by 45 right offset 0.5,0
set xtics nomirror
set grid xtics ytics lc rgb "#dddddd" lt 1
set terminal pngcairo font "sans,13" size 1920,1280
set bmargin 7;set lmargin 13;set rmargin 13;set tmargin 5
set ylabel "log(New {VV[1]} per day / New {VV[0]} per day)"

set output "{graphfn}"
set title "New cases per day in the UK of {VV[1]} compared with {VV[0]}\\nProgram: https://github.com/alex1770/Covid-19/blob/master/VOCgrowth/uk\\\\_var\\\\_comp.py\\nSource: COG-UK {cogdate}"
plot "{datafn}" u 1:(log($3/$2)):(sqrt($2*$3/($2+$3))/5) pt 5 ps variable title "log(Daily {VV[1]} / Daily {VV[0]}); larger blobs indicate more certainty (more samples)", (x/86400-{int(mindate)+cross})*{grad} lw 2 w lines title "{growthstr}\\n{doubstr}\\n{crossstr}"
"""

po=subprocess.Popen("gnuplot",shell=True,stdin=subprocess.PIPE)
p=po.stdin
p.write(cmd.encode('utf-8'))
p.close()
po.wait()
