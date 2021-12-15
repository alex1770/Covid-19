import sys
from math import log,exp
from scipy.stats import gamma as gammadist
from stuff import *

# Given al,be, need to specifiy three of g, V, NI, R, R0npi.
# Decided that the first three of these, g, V, NI, are the most reliable / directly observable, so write the sim function in terms of these.
# Results don't depend much on sd given g, mgt (i.e., spread in infectivity distribution not actually important given g).

def sim(date0="2021-07-12",
        mgt=5.5, sd=3,# parameters for Gamma infectivity distribution
        step4=1.5,# Multiplicative change in R0npi on 19 July (net effect of lifting of restrictions and school holidays)
        NI=16e6,# Number infected to date (maybe downrate those infected in first wave due to waning)
        g=0.04,
        newconfcases=30e3,car=0.4,# Current new confirmed cases per day, and estimated case ascertainment rate
        N=66.7e6,# Total population (ONS)
        dV=150e3,# Fully-vaccinate this number of people per day (i.e., dV = 1/2 daily dose rate)
        V0=40.2e6,# The number of people fully vaccinated (ONS)
        Vmax=53e6,# Vaccinatable population in UK (53m = 18+ ONS; 57.4m = 12+ ONS; both estimates)
        vhes=0.1,# Vaccine hesitancy rate in vaccinatable population
        veff=0.85,# Vaccine efficacy against transmission (includes efficacy against infection)
        ieff=0.8,# Infection efficacy against transmission (ditto)
        pr=False):
  
  # Get infectivity distribution from infectivity distribution Gamma(al,be)
  # id[i] = Expected number of infections on day i in unimmune population with unit R
  al=(mgt/sd)**2;be=al/mgt
  id=[0]
  i=1
  tp=0;te=0
  while 1:
    p=gammadist.cdf(be*(i+.5),al)-gammadist.cdf(be*(i-.5*(1+(i==1))),al)
    if i>al/be and p<1e-4: break
    id.append(p)
    i+=1
  nid=len(id)

  # R=(1+g/be)^al
  # Susceptible = N*(1-V/N)*(1-NI/N) = N*R/R0npi
  V=V0# Current number vaccinated
  S=N*(1-veff*V/N)*(1-ieff*NI/N)
  R=(1+g/be)**al
  R0npi=N*R/S# R-number with current NPIs, but with no immunity
  print("R0npi =",R0npi)
  I=[newconfcases*exp(g*i) for i in range(2-nid,1)]

  if pr: fp=sys.stdout
  else: fp=open('tempg','w')
  day0=datetoday(date0)
  TI=0
  while 1:
    n=len(I)
    day=day0+n-nid
    date=daytodate(day)
    if date=="2021-07-19": R0npi*=step4
    if date>="2021-09-10" and date<"2021-09-30": dV+=10e3# Assume new vaccine supply in September
    S=N*(1-veff*V/N)*(1-ieff*NI/N)
    t=0
    for i in range(1,nid): assert n-i>=0;t+=I[n-i]*id[i]*R0npi*S/N
    TI+=t
    V=min(V+dV,Vmax*(1-vhes))
    NI=min(NI+t,N)
    I.append(t)
    print("%s   %8.0f   %8.0f   %8.0f  %8.0f"%(date,t,NI,TI,Vmax*(1-vhes)-V),file=fp)
    if t<1 or date>="2022-03-01": break
  if not pr: fp.close()
  
  return TI

if 0:
  for vhes in [0.05, 0.1, 0.15]:
    for Vmax in [53e6,57e6]:
      for g in [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08]:
        print("g=%.2f   Vmax=%4.1fm  vhes=%4.2f   %4.1fm"%(g,Vmax/1e6,vhes,sim(g=g,Vmax=Vmax,vhes=vhes)/1e6))

print(sim(g=0.04,step4=1.3,Vmax=53e6,vhes=0.1,veff=0.85,ieff=0.8))

