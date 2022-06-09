#   Regard the generation time as a random variable, T.
#   Let J(t) denote the rate of new infections at time t (J=dI/dt), S(t)=susceptibles at t. Then:
#   J(t) = S(t)/N.R0.E[J(t-T)]
#   dS(t)/dt = -J(t)
#   
#   Get rid of the derivative:
#   dS(t)/dt = -J(t)
#            = -S(t)/N.R0.E[J(t-T)]
#            = S(t)/N.R0.E[dS(t-T)/dt]
#            = S(t)/N.R0.(d/dt)E[S(t-T)]
#   dlog(S(t))/dt = R0/N.(d/dt)E[S(t-T)]
#   log(S(t)) = R0/N.E[S(t-T)] + const
#
#   Here we look at the special case where gen time T=constant (=1 for convenience)
#   s(t)=S(t)/N, j(t)=J(t)/N
#   s(t) = C.exp(R0.s(t-1))
#   j(t) = -s'(t) = R0.s(t).j(t-1)
#   C=s(1)exp(-R0.s(0)) = s*.exp(-R0.s*), where s* is the upper susceptibility fixed point.
#   
#   s_{r+1} = C.exp(R0.s_r)
#   j_{r+1} = s_{r+1}.R0.j_r
#
#   You could choose s(t) = anything on [0,1) and then extend to [n,n+1), but I'm assuming we want
#   the real-analytic solution which is monotone decreasing. (Is that unique?)
#
#   Issues:
#   1) Ideally you'd parameterise by s(0), s'(0), but then calculating C isn't trivial (though it's easy to get a good approximation).
#   2) You'd like to know the intermediate values at every time step, not just at multiples of generation time (integers here).
#
#   These two issues are solved below, first to various approximations ("temp1" should be quick, and pretty accurate, apart from
#   extreme parameters, possibly), and finally properly/perfectly (with certain machine accuracy constraints).
#   The proper solution shifts time back until s is very close to the upper fixed point, s*, at which point
#   you can expand the iterated exponential around s* in this fashion:
#   https://en.wikipedia.org/wiki/Iterated_function#Some_formulas_for_fractional_iteration
#   s(t0+t) = s*-(s*-s(t0))*t^{R0.s*} + O((s*-s(t0))^2)
#   s'(t0)  = -(s*-s(t0))*log(R0.s*)
#   then can use iteration map to project the derivative s'(t0) forwards in time back to the present, where s=s0,
#   and compare the mapped derivative s'(t0) with desired s1.
#   This whole procedure was a function of C, so can choose C to match the derivative at s1.
#   If use n steps, and if f() is the inverse iteration map, implicitly depending on C and R0, then f(s)=log(s/C)/R0, and
#   the mapped derivative in the present is R0^n.s0.f(s0).f^2(s0).....f^{n-1}(s0).(s*-f^n(s0)).
#
#   Alternative to C: D=-log(C)/R0

from math import exp,log,sqrt

# Time is in units of generation time
s0=0.6# Initial number of susceptibles
s1=0.1# Initial new infection rate (rate of decline of susceptibles)
R0=1.1/s0

def F(C,s): return C*exp(R0*s)
def Fi(C,s): return log(s/C)/R0
def f(C,s): return C*R0*exp(R0*s)# dF(C,s)/ds
def G(D,s): return exp(R0*(s-D))
def Gi(D,s): return log(s)/R0+D
def g(D,s): return R0*exp(R0*(s-D))# dG(D,s)/ds

M=100# Subdivision of generation time
assert M%2==0

# SIR for comparison
s=s0;i=s1/(R0*s);r=0
gamma=1/M
tj=0
# d/dt(log(s)) = -R0.gamma.i
# d/dt(log(i)) = gamma.(R0.s-1)
with open('tempsir','w') as fp:
  while R0*s>1 or j>=1e-10:
    for it in range(M):
      j=R0*s*i
      print("%12.4f %10.6f %12g"%(r+it/M,s,j),file=fp)
      tj+=j/M
      s*=exp(-R0*gamma*i)
      i*=exp(gamma*(R0*s-1))
    r+=1
print(tj,s0-s)

# Renewal G, at multiples of generation time only, using estimated C, initial [-1/2,1/2]
C=(s0-s1/2)*exp(-R0*(s0+s1/2))# Approx
tj=0
with open('temp0','w') as fp:
  s=s0;j=s1;r=0
  while s*R0>1 or j>=1e-10:
    print("%4d %10.6f %12g"%(r,s,j),file=fp)
    s=F(C,s)
    tj+=j
    j*=s*R0
    r+=1
print(tj,s0-s)

# Renewal G, with intermediates, using estimated C, and a particular method for j (one of a few possibilities)
# s2 chosen to match j(1/2+)=j(1/2-) (and so all j(r+1/2+)=j(r+1/2-))
tj=0
X=R0*(s0-s1/2)
s2=(X-1)/(X+1)*2*s1
C=(s0-s1/2)*exp(-R0*(s0+s1/2))# Approx
with open('temp1','w') as fp:
  s_=[];j_=[]
  for r in range(M):
    t=r/M-0.5
    s_.append(s0-s1*t)
    j_.append(s1+s2*t)
  k=M//2;t=0
  while s_[0]*R0>1 or max(j_)>=1e-10:
    print("%12.4f %10.6f %12g"%(t,s_[k],j_[k]),file=fp)
    tj+=j_[k]/M
    k+=1;t+=1/M
    if k==M:
      s_=[F(C,s) for s in s_]
      j_=[j*s*R0 for j,s in zip(j_,s_)]
      k=0
print(tj,s0-s_[0])

# Renewal G, with intermediates to order 2, using estimated C, initial [0,1]
# Should have used [-1/2,1/2] which is totally doable but a big PITA.
tj=0
QA=R0/2
QB=R0*s0+(1/2)*R0**2*s1**2+1
QC=((R0*s1+2)*R0*(s0-s1)-2)*s1
s2=(QB-sqrt(QB**2-4*QA*QC))/(2*QA)
c0=s0-s1-s2/2
C=c0*exp(-R0*s0)# Approx
s3=2*((R0*c0-1)*s1-s2)
with open('temp2','w') as fp:
  s_=[];j_=[]
  for r in range(M):
    t=(r+0.5)/M
    s_.append(s0-s1*t-s2*t**2/2)
    j_.append(s1+s2*t+s3*t**2/2)
  k=0
  while s_[0]*R0>1 or max(j_)>=1e-10:
    for r in range(M):
      print("%12.4f %10.6f %12g"%(k+(r+0.5)/M,s_[r],j_[r]),file=fp)
    tj+=sum(j_)/M
    s_=[F(C,s) for s in s_]
    j_=[j*s*R0 for j,s in zip(j_,s_)]
    k+=1
print(tj,s0-s)

# Renewal G, exact (in limit where take small epilsons to 0), though may need to robustify
# Using D=-log(C)/R0 instead of C, because if s1 is large then C can get too small.

def upperfixedpoint(D):
  st=1
  while 1:
    prev=st
    st=Gi(D,st)
    if abs(st-prev)<1e-12: break
  return st  

# Go back to early time when everything is linear, calculate the derivative s'(t) there, and the project forward to the present
def sl(D):
  st=upperfixedpoint(D)
  pr=log(R0*st);s=s0
  while 1:
    pr*=R0*s
    s=Gi(D,s)
    #print(s,st-s,pr*(st-s))
    if st-s<1e-5: break
  return pr*(st-s)

minD=(log(R0)+1-1e-8)/R0

# First approximation to D
D=(-log(s0)+s1/s0/2+R0*s0*exp(s1/s0/2))/R0
D=max(D,minD)
# D=minD

# Find some values of D bracketing the desired one (though it's not actually necessary to bracket)
D0=D1=D
sl0=sl1=sl(D)
while (sl1<s1)==(sl0<s1):
  if sl1<s1: D1+=0.1
  else: D1=max(D1-0.1,minD)
  sl1=sl(D1)

# Find D solving sl(D)=s1 by a Newton-type method
while 1:
  D=(s1-sl0)/(sl1-sl0)*(D1-D0)+D0
  if abs(D1-D0)<1e-6: break
  sl2=sl(D)
  if abs(sl0-s1)<abs(sl1-s1): D1=D;sl1=sl2
  else: D0=D;sl0=sl2

# Go back to early time near the fixed point
st=upperfixedpoint(D)
s=s0;it=0
while 1:
  s=Gi(D,s)
  it+=1
  if st-s<1e-5: break

# Use iterated exponential expansion near the fixed point
# lam is chosen such that
# G(D,st-(st-s)*exp(-lam/2))=st-(st-s)*exp(lam/2)
# Expect lam~=log(R0*st)
lam=log(R0*st)
err=G(D,st-(st-s)*exp(-lam/2))-(st-(st-s)*exp(lam/2))
err1=g(D,st-(st-s)*exp(-lam/2))*(st-s)*exp(-lam/2)/2-(-(st-s)*exp(lam/2)/2)
lam-=err/err1
s_=[]
for k in range(M+1):
  t=k/M-0.5
  #s_.append(st-(st-s)*(R0*st)**t)
  s_.append(st-(st-s)*exp(lam*t))

# Project forward to the present
for j in range(it):
  s_=[G(D,s) for s in s_]

# Print out results
srecent=s_[M//2-1],s_[M//2],s_[M//2+1]
t=0
k=M//2+1
tj=0
with open('temp3','w') as fp:
  while 1:
    s=srecent[1]
    j=-(srecent[2]-srecent[0])/2*M
    tj+=j/M
    print("%12.4f %10.6f %12g"%(t,s,j),file=fp)
    k+=1
    if k==M:
      k=0
      s_=[G(D,s) for s in s_]
    s=s_[k]
    srecent=(srecent[1],srecent[2],s)
    t+=1/M
    if not (s*R0>1 or j>=1e-10): break
print(tj,s0-s)
