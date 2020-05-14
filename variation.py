# Exploring https://www.medrxiv.org/content/10.1101/2020.04.27.20081893v1

from scipy.special import gammainc
from scipy.stats import gamma as gammadist
import numpy as np
from math import sqrt

# Obtained (this version of) the number of initial infections by adding the number of new
# cases (as reported by Worldometers) on the days Feb 27, 28, 29, Mar 1, then multiplying
# by 10. Four days, because 1/gamma=4, and multiplying by 10 because p=0.1 (assumed
# probability of an infection being reported as a case).

situations=[
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "gamma",
    "CV": 1,
    "SD": 0.00
  },
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "gamma",
    "CV": 3,
    "SD": 0.00
  },
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "gamma",
    "CV": 1,
    "SD": 0.68
  },
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "gamma",
    "CV": 3,
    "SD": 0.61
  },
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "twopoint",
    "CV": 1,
    "x": 0,
    "SD": 0.68
  },
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "twopoint",
    "CV": 1,
    "x": 0.99,
    "SD": 0.68
  },
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "twopoint",
    "CV": 3,
    "x": 0,
    "SD": 0.61
  },
  {
    "Country": "Italy",
    "Population": 60360000,
    "initialinfections": 12310,
    "disttype": "twopoint",
    "CV": 3,
    "x": 0.99,
    "SD": 0.61
  },
  {
    "Country": "Austria",
    "Population": 8859000,
    "initialinfections": 120,
    "disttype": "gamma",
    "CV": 1,
    "SD": 0.00
  },
  {
    "Country": "Austria",
    "Population": 8859000,
    "initialinfections": 120,
    "disttype": "gamma",
    "CV": 3,
    "SD": 0.00
  },
  {
    "Country": "Austria",
    "Population": 8859000,
    "initialinfections": 120,
    "disttype": "gamma",
    "CV": 1,
    "SD": 0.80
  },
  {
    "Country": "Austria",
    "Population": 8859000,
    "initialinfections": 120,
    "disttype": "gamma",
    "CV": 3,
    "SD": 0.77
  }
]

delta=1/4
gamma=1/4
rho=0.5
R0=2.7
p=0.1
days=487# Days from 2020-03-01 to 2021-07-01

# By experimentation, values of 10 for stepsperday and sbins give near to limiting
# behaviour (within small tolerances), so 100 for each is hopefully plenty.
stepsperday=100# Subdivision of each day
maxsbins=100# Number of bins for susceptibility (equally spaced by CDF)

# Caption to Figure 1 gives this time-dependence of the social distancing parameter:
# Return value, x, is effective social distancing number. Infection force is multiplied
# by 1-x.
def SD(sd0,day):
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

for sit in situations:
  N=sit["Population"]
  cv=sit["CV"]
  sd0=sit["SD"]
  print("Country:",sit["Country"])
  print("Coefficient of Variation:",cv)
  print("Max social distancing:",sd0)

  susc,q,desc = getsusceptibilitydist(sit,maxsbins)
  print("Susceptibility distribution:",desc)
  sbins=len(susc)
  checkcv(susc,q,cv)
  S=q*N

  if cv not in seen:
    seen.add(cv)
    with open('q_CV%g'%cv,'w') as fp:
      for i in range(sbins):
        print("%9.3f   %9.7f"%(susc[i],q[i]),file=fp)
  
  E=np.zeros(sbins)
  I=np.zeros(sbins)
  beta=R0/(rho/delta+1/gamma)# NB there is a factor of 1/N error in formula (2) from the paper
  
  lam=0

  fn='output_%s_%s_SD%g'%(sit['Country'],desc,sd0)
  with open(fn,'w') as fp:
    print("#   Day            s        e        i   I_reported",file=fp)
    for d0 in range(days*stepsperday):
      day=d0/stepsperday
  
      # Introduce initial infections by boosting the infection force for the expected duration of an infection (1/gamma days)
      if day<1/gamma: lam+=beta/N*sit['initialinfections']
      
      new=lam*susc*S*(1-SD(sd0,day))
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

