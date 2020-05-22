# Exploring https://www.medrxiv.org/content/10.1101/2020.04.27.20081893v1
# This reproduces the susceptibilty curves for v1 of the paper, but not the connectivity curves (possibly a discrepancy in initial data)

from scipy.special import gammainc
from scipy.stats import gamma as gammadist
from scipy.stats import norm as normaldist
import numpy as np
from math import sqrt,log

# Set initial infections to replicate the output in the Gomes paper (checking they are plausible values).

countries=[
  {
    "name": "Italy",
    "population": 60360000,
    "initialinfections": 14000,
    "delay": 10# Number of days before phasing in lockdown
  },
  {
    "name": "Austria",
    "population": 8859000,
    "initialinfections": 540,
    "delay": 16# Number of days before phasing in lockdown
  }
]

distributions=[
  {
    "disttype": "gamma",
    "CV": 1,
  },
  {
    "disttype": "gamma",
    "CV": 3,
  },
  {
    "disttype": "lognormal",
    "CV": 1,
  },
  {
    "disttype": "lognormal",
    "CV": 3,
  },
  {
    "disttype": "twopoint",
    "CV": 1,
    "x": 0,
  },
  {
    "disttype": "twopoint",
    "CV": 1,
    "x": 0.9,
  },
  {
    "disttype": "twopoint",
    "CV": 3,
    "x": 0,
  },
  {
    "disttype": "twopoint",
    "CV": 3,
    "x": 0.9,
  }
]

SD={
  "Italy": {
    "susceptibility": {1: 0.68, 3: 0.61},
    "connectivity": {1: 0.67, 3: 0.65}
  },
  "Austria": {
    "susceptibility": {1: 0.80, 3: 0.77},
    "connectivity": {1: 0.80, 3: 0.78}
  }
}

# R0 values used in v1 of the paper, as a function of mode (susceptibility/connectivity) and CV:
R0values={
  "susceptibility": {1: 2.7, 3: 2.7},
  "connectivity": {1: 2.8, 3: 3.1}
}

delta=1/4# rate of progression (in days^{-1}) from E->I
gamma=1/4# rate of progression (in days^{-1}) from I->R
rho=0.5#   Relative infectivity of E group compared with I group
p=0.026#   Proportion of infections that are reported. The paper says to use 0.1, but 0.026 seems to replicate their
#          results much better in the susceptibility mode.
days=487#  Days from 2020-03-01 to 2021-07-01

# By experimentation, values of 10 for stepsperday and 50 for sbins give near to limiting
# behaviour (within small tolerances), so 100 and 500 is hopefully plenty.
stepsperday=100# Subdivision of each day
maxsbins=500# Number of bins for susceptibility (equally spaced by CDF)

# Caption to Figure 1 give this time-dependence of the social distancing parameter:
# return value, x, is effective social distancing number (infection force is multiplied by 1-x).
# Initial 10 or 16 day period ("delay") is from examining the red R0 and Rt graphs.
def SDt(day,delay):
  if day<delay: return 0
  day-=delay
  if day<14: return day/14
  day-=14
  if day<31: return 1
  day-=31
  if day<365: return 1-day/365
  return 0

# Returns template distribution with mean 1 and coefficient of variation cv, using up to sbins bins
def getsusceptibilitydist(sit,maxsbins):

  cv=sit["CV"]
  
  if sit["disttype"]=="gamma":
    # Shape parameter of Gamma distribution
    k=1/cv**2
    
    # Space out bins according to shape k+1 so that each bin represents an equal
    # (incomplete) expectation of X in the shape k distribution. Other bin choices are
    # possible, but this is a good one from the point of view of getting a more accurate
    # CV using a smaller number of bins.
    l=gammadist.ppf([i/maxsbins for i in range(maxsbins)],k+1)
  
    # Calculate:
    #   m0[i]   = P[X < l_i]
    #   m1[i]   = E[X; X < l_i]
    #   susc[i] = E[X | l_i <= X < l_{i+1}], the representative susceptibility for bin i
    #   q[i]    = P(l_i <= X < l_{i+1})
    m0=np.append(gammainc(k,l),1)
    m1=np.append(gammainc(k+1,l),1)
    susc=np.array([(m1[i+1]-m1[i])/(m0[i+1]-m0[i]) for i in range(maxsbins)])
    q=np.array([m0[i+1]-m0[i] for i in range(maxsbins)])
    desc="gamma_CV%g"%cv

  elif sit["disttype"]=="twopoint":
    x=sit["x"]
    p=cv**2/((1-x)**2+cv**2)
    y=(1-p*x)/(1-p)
    q=np.array([p,1-p])
    susc=np.array([x,y])
    desc="twopoint_CV%g_x%g"%(cv,x)

  elif sit["disttype"]=="lognormal":
    var=log(1+cv**2)
    mu=-var/2

    # Space out bins according to lognormal(-3mu,var) so that each bin represents an equal
    # (incomplete) expectation of X^2 in the lognormal(mu,var) distribution. Other bin
    # choices are possible, but this is a good one from the point of view of getting a
    # more accurate CV using a smaller number of bins.
    l=normaldist.ppf([i/maxsbins+1e-30 for i in range(maxsbins)],-3*mu,sqrt(var))

    # Given
    #   L_i    = exp(l_i)
    #   mu     = -var/2
    #   log(X) ~ N(mu,var)
    # Calculate:
    #   m0[i]   = P[X < L_i] = P[log(X) < l_i]
    #   m1[i]   = E[X; X < L_i] = E[X; log(X) < l_i]
    #   susc[i] = E[X | L_i <= X < L_{i+1}], the representative susceptibility for bin i
    #   q[i]    = P(L_i <= X < L_{i+1})
    m0=np.append(normaldist.cdf(l,mu,sqrt(var)),1)
    m1=np.append(normaldist.cdf(l,-mu,sqrt(var)),1)
    susc=np.array([(m1[i+1]-m1[i])/(m0[i+1]-m0[i]) for i in range(maxsbins)])
    q=np.array([m0[i+1]-m0[i] for i in range(maxsbins)])
    desc="lognormal_CV%g"%cv

  else: raise NotImplementedError("Unknown susceptibility distribution type: "+sit["disttype"])
  
  return susc,q,desc
  

def checkcv(susc,q,cv):
  x0=x1=x2=0
  err=0
  for (v,p) in zip(susc,q):
    x0+=p
    x1+=p*v
    x2+=p*v**2
  mu=x1/x0
  sd=sqrt(x2-x1**2)
  err=max(abs(x0-1),abs(mu-1),abs(sd/mu-cv))
  print("Distribution quantisation error: %g"%err)

seen=set()

for country in countries:
  N=country['population']
  for mode in ['susceptibility', 'connectivity']:
    for dist in distributions:
      cv=dist["CV"]
      R0=R0values[mode][cv]
      for sd0 in [0, SD[country['name']][mode][cv]]:
        #if not (country['name']=="Italy" and mode=="connectivity" and cv==3 and sd0>=0.64): continue
        print("Country:",country['name'])
        print("Mode:",mode)
        print("R0:",R0)
        print("Coefficient of Variation:",cv)
        print("Max social distancing:",sd0)
      
        susc,q,desc = getsusceptibilitydist(dist,maxsbins)
        print("Susceptibility distribution:",desc)
        sbins=len(susc)
        checkcv(susc,q,cv)
        S=q*N
      
        if desc not in seen:
          seen.add(desc)
          with open('q_%s'%desc,'w') as fp:
            for i in range(sbins):
              print("%9.3f   %9.7f"%(susc[i],q[i]),file=fp)

        # NB there is a factor of 1/N error in formula (2) from the paper
        if mode=='susceptibility':
          beta=R0/(rho/delta+1/gamma)
        else:
          beta=R0/((1+cv**2)*(rho/delta+1/gamma))
        E=np.zeros(sbins)
        I=np.zeros(sbins)
        # Assume initial infections occur proportional to susceptibility by including the factor susc[i] here
        # (though this won't make a big difference)
        for i in range(sbins):
          E[i]=country['initialinfections']*q[i]*susc[i]
          S[i]-=I[i];assert S[i]>=0
        
        fn='output_%s_%s_%s_SD%g'%(country['name'],mode,desc,sd0)
        HIT=None
        with open(fn,'w') as fp:
          print("#   Day            s        e        i   I_reported      R_s     R_t",file=fp)
          for d0 in range(days*stepsperday):
            day=d0/stepsperday
            Ssum=S.sum()
            Esum=E.sum()
            Isum=I.sum()
            if mode=='connectivity':
              Esum_infective=(susc*E).sum()
              Isum_infective=(susc*I).sum()
              R_s=(susc*susc*S).sum()*beta/N*(rho/delta+1/gamma)
            else:
              Esum_infective=Esum
              Isum_infective=Isum
              R_s=(susc*S).sum()*beta/N*(rho/delta+1/gamma)
            sd=sd0*SDt(day,country['delay'])# current social distancing
            R_t=R_s*(1-sd)
            if HIT==None and R_s<=1: HIT=1-Ssum/N
            print("%7.2f      %7.5f  %7.5f  %7.5f    %9.0f   %6.3f  %6.3f"%(day,Ssum/N,Esum/N,Isum/N,p*Isum,R_s,R_t),file=fp)
            lam=beta/N*(rho*Esum_infective+Isum_infective)
            new=lam*susc*S*(1-sd)
            I+=(delta*E-gamma*I)/stepsperday
            E+=(new-delta*E)/stepsperday
            S+=-new/stepsperday
          final=(1-Ssum/N)*100
          print("Herd immunity threshold ",end="")
          if HIT==None: print("> %.1f%% (not yet attained)"%final)
          else: print("= %.1f%%"%(HIT*100))
          print("Final proportion infected ",end="")
          if lam>1e-4: print("> %.1f%% (infection ongoing)"%final)
          else: print("= %.1f%%"%final)
        print("Written output to file \"%s\""%fn)
        print()

# gnuplot> plot "output_Italy_susceptibility_gamma_CV1_SD0" u 1:5 w lines, "output_Italy_susceptibility_gamma_CV1_SD0.68" u 1:5 w lines
# gnuplot> plot "output_Italy_susceptibility_gamma_CV3_SD0" u 1:5 w lines, "output_Italy_susceptibility_gamma_CV3_SD0.61" u 1:5 w lines
# gnuplot> plot "output_Italy_connectivity_gamma_CV1_SD0" u 1:5 w lines, "output_Italy_connectivity_gamma_CV1_SD0.67" u 1:5 w lines
# gnuplot> plot "output_Italy_connectivity_gamma_CV3_SD0" u 1:5 w lines, "output_Italy_connectivity_gamma_CV3_SD0.65" u 1:5 w lines 
# etc.
