# Exploring https://www.medrxiv.org/content/10.1101/2020.04.27.20081893v1

from scipy.special import gammainc
from scipy.stats import gamma as gammadist
import numpy as np
from math import sqrt

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
    "disttype": "twopoint",
    "CV": 1,
    "x": 0,
  },
  {
    "disttype": "twopoint",
    "CV": 1,
    "x": 0.99,
  },
  {
    "disttype": "twopoint",
    "CV": 3,
    "x": 0,
  },
  {
    "disttype": "twopoint",
    "CV": 3,
    "x": 0.99,
  }
]

SD={
  "Italy": {1: 0.68, 3: 0.61},
  "Austria": {1: 0.80, 3: 0.77}
}

delta=1/4# rate of progression (in days^{-1}) from E->I
gamma=1/4# rate of progression (in days^{-1}) from I->R
rho=0.5# Relative infectivity of E group compared with I group
R0=2.7
p=0.026# Proportion of infections that are reported. The paper says to use 0.1, but 0.026 seems to replicate their
#        results much better.
days=487# Days from 2020-03-01 to 2021-07-01

# By experimentation, values of 10 for stepsperday and sbins give near to limiting
# behaviour (within small tolerances), so 100 for each is hopefully plenty.
stepsperday=100# Subdivision of each day
maxsbins=100# Number of bins for susceptibility (equally spaced by CDF)

# Caption to Figure 1 give this time-dependence of the social distancing parameter:
# return value, x, is effective social distancing number (infection force is multiplied by 1-x).
# Initial 10 or 16 day period ("delay") is from examining the red R0 and Rt graphs.
def SDt(sd0,day,delay):
  if day<delay: return 0
  day-=delay
  if day<14: return day/14*sd0
  day-=14
  if day<31: return sd0
  day-=31
  if day<365: return (1-day/365)*sd0
  return 0

# Returns template distribution with mean 1 and coefficient of variation cv, using up to sbins bins
def getsusceptibilitydist(sit,maxsbins):

  cv=sit["CV"]
  
  if sit["disttype"]=="gamma":
    # Shape parameter of Gamma distribution
    k=1/cv**2
    
    # Space out bins according to shape k+1 so that each bin represents an equal
    # (incomplete) expectation. Other bin choices are possible, but k+1 is better than using
    # k (which would make each bin represent an equal probability from the gamma
    # distribution) in the sense that you get limiting behaviour for a smaller value of
    # maxsbins. (This is to be expected because the infection force is the important quantity.)
    l=gammadist.ppf([i/maxsbins for i in range(maxsbins)],k+1)
  
    # Calculate:
    #   m0[i] = P[X<i/sbin]
    #   m1[i] = E[X; X<i/sbin]
    #   susc[i] = E[X | i/sbin < X < (i+1)/sbin], the representative susceptibility for bin i
    #   q[i] = P(i/sbin < X < (i+1)/sbin)
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

  else: raise NotImplementedError("Unknown susceptibility distribution type "+sit["disttype"])
  
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
  for dist in distributions:
    cv=dist["CV"]
    for sd0 in [0, SD[country['name']][cv]]:
      print("Country:",country['name'])
      print("Coefficient of Variation:",cv)
      print("Max social distancing:",sd0)
    
      susc,q,desc = getsusceptibilitydist(dist,maxsbins)
      print("Susceptibility distribution:",desc)
      sbins=len(susc)
      checkcv(susc,q,cv)
      S=q*N
    
      if cv not in seen:
        seen.add(cv)
        with open('q_CV%g'%cv,'w') as fp:
          for i in range(sbins):
            print("%9.3f   %9.7f"%(susc[i],q[i]),file=fp)
      
      beta=R0/(rho/delta+1/gamma)# NB there is a factor of 1/N error in formula (2) from the paper
      E=np.zeros(sbins)
      I=np.zeros(sbins)
      # Assume initial infections occur proportional to susceptibility by including the factor susc[i] here
      # (though this won't make a big difference)
      for i in range(sbins):
        I[i]=country['initialinfections']*q[i]*susc[i]
        S[i]-=I[i];assert S[i]>=0
      
      lam=0
    
      fn='output_%s_%s_SD%g'%(country['name'],desc,sd0)
      with open(fn,'w') as fp:
        print("#   Day            s        e        i   I_reported",file=fp)
        for d0 in range(days*stepsperday):
          day=d0/stepsperday
      
          new=lam*susc*S*(1-SDt(sd0,day,country['delay']))
          I+=(delta*E-gamma*I)/stepsperday
          E+=(new-delta*E)/stepsperday
          S+=-new/stepsperday
          Ssum=S.sum()
          Esum=E.sum()
          Isum=I.sum()
          lam=beta/N*(rho*Esum+Isum)
          print("%7.2f      %7.5f  %7.5f  %7.5f    %9.0f"%(day,Ssum/N,Esum/N,Isum/N,p*Isum),file=fp)
        print("Final proportion infected = %.1f%%"%((1-Ssum/N)*100))
      print("Written output to file \"%s\""%fn)
      print()

# gnuplot> plot "output_Italy_gamma_CV1_SD0" u 1:5 w lines, "output_Italy_gamma_CV1_SD0.68" u 1:5 w lines
# gnuplot> plot "output_Italy_gamma_CV3_SD0" u 1:5 w lines, "output_Italy_gamma_CV3_SD0.61" u 1:5 w lines
# etc.
