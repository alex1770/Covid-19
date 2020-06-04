
from __future__ import division, print_function

# Evaluate the predictive power of poo on hospital admissions
# Re-evaluating this paper https://www.medrxiv.org/content/10.1101/2020.05.19.20105999v1
# following Nick Brown's critique, and using his extracted data
# https://twitter.com/sTeamTraen/status/1265411882283917315
# Data from https://twitter.com/sTeamTraehttps://twitter.com/sTeamTraen/status/1265411887862349826/photo/1n/status/1265411887862349826/photo/1

poo=[5,5,25,0,50,30,80,80,85,75,90,95,105,45,110,70,140,115,175,105,130,125,290,120,350,185,135,205,195,80,65,100,90,120,45,55,45,90,80,85,25,55,85,90]
hosp=[8,10,3,8,10,12,8,10,17,14,14,26,17,12,27,21,15,17,15,21,22,21,22,25,10,29,27,20,26,25,20,5,14,12,18,23,12,17,13,15,18,11,19,19]

n=len(poo)
assert len(hosp)==n

from math import sqrt

# Return correlation coefficient
def corr(xx,yy):
  n=len(xx);assert len(yy)==n
  mux=sum(xx)/n;muy=sum(yy)/n
  sxy=sum((x-mux)*(y-muy) for (x,y) in zip(xx,yy))
  sxx=sum((x-mux)*(x-muy) for x in xx)
  syy=sum((y-mux)*(y-muy) for y in yy)
  return sxy/sqrt(sxx*syy)

# Return adjusted R^2 for xx ==> yy
def R2_1(X,Y):
  n=len(X);assert len(Y)==n
  mux=sum(X)/n;muy=sum(Y)/n
  xx=xy=yy=0
  for (x,y) in zip(X,Y):
    x-=mux;y-=muy
    xx+=x*x
    xy+=x*y
    yy+=y*y
  al=xy/xx
  resid=al**2*xx-2*al*xy+yy
  return 1-(resid/(n-2))/(yy/(n-1))

# Return adjusted R^2 for xx,yy ==> zz
def R2_2(X,Y,Z):
  n=len(X);assert len(Y)==n and len(Z)==n
  mux=sum(X)/n;muy=sum(Y)/n;muz=sum(Z)/n
  xx=xy=yy=xz=yz=zz=0
  for (x,y,z) in zip(X,Y,Z):
    x-=mux;y-=muy;z-=muz
    xx+=x*x
    xy+=x*y
    yy+=y*y
    xz+=x*z
    yz+=y*z
    zz+=z*z
  det=xx*yy-xy**2
  al=(yy*xz-xy*yz)/det
  be=(-xy*xz+xx*yz)/det
  resid=al**2*xx+be**2*yy+zz+2*al*be*xy-2*al*xz-2*be*yz
  return 1-(resid/(n-3))/(zz/(n-1)), al, be

print("Simple time-shifted correlations (h_i with h_{i+d} and p_i with h_{i+d})")
for day in range(-10,11):
  # poo[i] vs hosp[i+day]
  i0=max(0,-day)
  i1=min(n,n-day)
  print("Day %+3d   corr_hh = %+5.3f   corr_ph = %+5.3f"%(day,corr(hosp[i0:i1],hosp[i0+day:i1+day]),corr(poo[i0:i1],hosp[i0+day:i1+day])))
print()
  
print("Compare adjusted R^2s (h_i predicting h_{i+d} cf h_i, p_i predicting h_{i+d})")
for day in range(-10,11):
  # poo[i] vs hosp[i+day]
  i0=max(0,-day)
  i1=min(n,n-day)
  R2_hh=R2_1(hosp[i0:i1],hosp[i0+day:i1+day])
  R2_phh,al,be=R2_2(poo[i0:i1],hosp[i0:i1],hosp[i0+day:i1+day])
  print("Day %+3d  R^2(h==>h') = %+5.3f   R^2(p,h==>h') = %+5.3f   Poo advantage = %+5.3f    coeff_poo = %+5.3f  coeff_hosp = %+5.3f"%(day,R2_hh,R2_phh,R2_phh-R2_hh,al,be))
print()
