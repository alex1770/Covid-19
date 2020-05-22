# Terribly inefficient, but just a one-off

from math import exp,log,sqrt

from scipy.stats import gamma, norm, lognorm

# CV -> (r,Lambda) -> E[X^r.exp(-Lambda.X)], X ~ Gamma(shape=1/CV^2,scale=CV^2)
def makegammaexpectation(CV):
  if CV>0:
    k=1/CV**2# shape
    s=CV**2# scale
    def gammaexpectation(r,Lambda):
      return gamma.expect(lambda x:x**r*exp(-Lambda*x),(k,),scale=s)
    return gammaexpectation
  else:
    def constexpectation(r,Lambda):
      return exp(-Lambda)
    return constexpectation

# CV -> (r,Lambda) -> E[X^r.exp(-Lambda.X)], X ~ exp(N(var,-var/2)), var=log(1+CV^2)
def makelognormalexpectation(CV):
  if CV>0: # lognormal can't cope with variance 0 for some reason
    var=log(1+CV**2)
    mu=-var/2
    def lognormalexpectation(r,Lambda):
      return lognorm.expect(lambda x:x**r*exp(-Lambda*x),(sqrt(var),),scale=exp(mu))
    return lognormalexpectation
  else:
    def constexpectation(r,Lambda):
      return exp(-Lambda)
    return constexpectation

# CV,x -> (r,Lambda) -> E[X^r.exp(-Lambda.X)], X ~ twopoint(CV,x)
def maketwopointexpectation(CV,x):
  p=CV**2/((1-x)**2+CV**2)
  y=(1-p*x)/(1-p)
  def twopointexpectation(r,Lambda):
    return p*x**r*exp(-Lambda*x)+(1-p)*y**r*exp(-Lambda*y)
  return twopointexpectation

def HIT(R0,dist,connectivity=False):
  e0=dist(0,0)
  e1=dist(1,0)
  if abs(e0-1)+abs(e1-1)>1e-6: raise ArithmeticError("Distribution not normalised with mean 1")
  l0=0
  if connectivity:
    l1=2*sqrt(R0)/exp(1)
    R=R0/dist(2,0)
    r=2
  else:
    l1=R0/exp(1)
    R=R0
    r=1
  while l1-l0>1e-6:
    l=(l0+l1)/2
    v=R*dist(r,l)
    if v>1: l0=l
    else: l1=l
  return 1-dist(0,(l0+l1)/2)

R0=3
x0=0
x1=0.9

fn="hitout"
with open(fn,"w") as fp:
  print("#    R0 =",R0,file=fp)
  print("#    CV     Gamma-susc     Gamma-conn   Lognorm-susc   Lognorm-conn      TP-0-susc      TP-0-conn    TP-0.9-susc    TP-0.9-conn",file=fp)
  for cv0 in range(401):
    CV=cv0/100
    print("%7.4f"%CV,end="",file=fp)
    l=[]
    for dist in (makegammaexpectation(CV), makelognormalexpectation(CV), maketwopointexpectation(CV,x0), maketwopointexpectation(CV,x1)):
      for conn in [False, True]:
        print("        %7.4f"%HIT(R0,dist,conn),end="",file=fp)
    print(file=fp)
    fp.flush()

print('Written output to file "'+fn+'"')
