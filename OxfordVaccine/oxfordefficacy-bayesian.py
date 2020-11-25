n=131
(p0,v0,p1,v1)=(30,3,71,27)

from scipy.stats import gamma

for (pl,va) in [(p0,v0), (p1,v1), (p0+p1,v0+v1)]:

  # Gamma prior in shape, rate notation
  f=1
  alpha=1/f
  beta=1/((pl+va)/2)*f
  
  # Posterior distribution is (independently)
  # lamba ~ Gamma(pl+alpha, 1+beta)
  # mu    ~ Gamma(va+alpha, 1+beta)
  # We're interested in 95% (say) point of the mu/lambda distribution
  
  N=1000000
  Lam=gamma.rvs(pl+alpha, scale=1/(1+beta), size=N)
  Mu= gamma.rvs(va+alpha, scale=1/(1+beta), size=N)
  l=Mu/Lam
  l.sort()

  print("Number of infected participants on placebo arm:",pl)
  print("Number of infected participants on vaccine arm:",va)
  print("Using prior of Gamma(shape=%f, scale=%f)"%(alpha,1/beta))
  for c in [0.1, 0.25, 0.5, 0.9, 0.95, 0.98, 0.99]:
    print("%4.0f%% credible that efficacy >= %4.1f%%"%(c*100,(1-l[int(c*N)])*100))
  print()

