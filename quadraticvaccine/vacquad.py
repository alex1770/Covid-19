from math import log,exp,sqrt
import numpy as np

# Two population (V=to-be-vaccinated, N=not-to-be-vaccinated) SIR with vaccine term (alpha)
# Assume a simple contact matrix C = k.[p  (1-p)(1-v)]
#                                      [(1-p)v      p]
#
# where v is the proportion of people in the vaccinated group (assumed small), 1/2 <= p <= 1, and k is a scale factor, chosen so that
# C has its largest eigenvalue equal to 1, so that R_0=beta/gamma as usual.
#
# Note that this is a definite assumption. The general normalised contact matrix would have three parameters
# (just restricting biggest eigenvalue to be 1) whereas we're only using two (p and v).
# It's because we're restricting to contact matrices where a member of V tends to meet
# a similar (for non-extreme p) number of members of V to those of N.
# (Which implies a member of N tends to meet only a small number of members of V, because v is small.)
#
#
# I_V, S_V, I_N, S_N represent proportions of their respective populations (so between 0 and 1)
#
# Force term:
#        [ F_V ] = beta.C.[ I_V ]
#        [ F_N ]          [ I_N ]
#
#
# (d/dt) [ I_V ] =   [ F_V.S_V ] - gamma.[ I_V ]
#        [ I_N ]     [ F_N.S_N ]         [ I_N ]
#  
# (d/dt) [ S_V ] = - [ F_V.S_V ] - alpha.[ 1 ]
#        [ S_N ]     [ F_N.S_N ]         [ 0 ]
#
# Aim to work out I_V, I_N to second order in t in order to get expressions of the form
# I_V(t) = I_V0*exp( L_V1*t + (1/2)L_V2*t^2 + O(t^3) )
# I_N(t) = I_N0*exp( L_N1*t + (1/2)L_N2*t^2 + O(t^3) )

v=3/67# Vaccinated group size as a proportion of the population (using 80+ as example: 3m/67m)
I_N0=0.02
I_V0=0.01
S_N0=0.8
S_V0=0.9
beta=0.09
gamma=0.09
alpha=0.01

print("#    p           L_V1          L_V2          L_N1          L_N2")

for p in [0.01*x for x in range(50,101)]:

  C=np.array([[p,(1-p)*(1-v)],[(1-p)*v,p]])
  e=np.sort(np.abs(np.linalg.eigvals(C)))[1]# Largest eigenvalue of C
  C/=e
  
  # First work out S_V, S_N to first order in t
  # S_V = S_V0 + S_V1.t + ...
  # S_N = S_N0 + S_N1.t + ...
  F0=beta*np.matmul(C,[I_V0,I_N0])
  [S_V1, S_N1] = -([S_V0,S_N0]*F0) - np.array([alpha,0])
  
  # (d/dt) [I_V] = ( beta.[ S_V     ].C - gamma.[ 1    ] ).[I_V] = (A + Bt + O(t^2)).[I_V]
  #        [I_N]          [     S_N ]           [    1 ]   [I_N]                     [I_N]
  #
  A=beta*np.matmul([[S_V0,0],[0,S_N0]],C) - gamma*np.array([[1,0],[0,1]])
  B=beta*np.matmul([[S_V1,0],[0,S_N1]],C)

  # Solving above DE gives:
  # [I_V] = (I + At + (1/2)(A^2+B)t^2 + ...) [I_V0]
  # [I_N]                                    [I_N0]
  # 
  # I_V = I_V0 + I_V1.t + (1/2)I_V2.t^2 + ...
  # I_N = I_N0 + I_N1.t + (1/2)I_N2.t^2 + ...
  #
  [I_V1, I_N1] = np.matmul(A,[I_V0,I_N0])
  [I_V2, I_N2] = np.matmul(np.matmul(A,A)+B,[I_V0,I_N0])

  # log(I_V) = log(I_V0) + I_V1/I_V0.t + (1/2)(I_V2/I_V0-(I_V1/I_V0)^2).t^2 + ...
  # log(I_N) = log(I_N0) + I_N1/I_N0.t + (1/2)(I_N2/I_N0-(I_N1/I_N0)^2).t^2 + ...
  #
  L_V0=log(I_V0);L_V1=I_V1/I_V0;L_V2=I_V2/I_V0-(I_V1/I_V0)**2
  L_N0=log(I_N0);L_N1=I_N1/I_N0;L_N2=I_N2/I_N0-(I_N1/I_N0)**2
  print("%6.3f   %12.8f  %12.8f  %12.8f  %12.8f"%(p,L_V1,L_V2,L_N1,L_N2))

  # Check
  if 0:
    (I_V,S_V,I_N,S_N)=(I_V0,S_V0,I_N0,S_N0)
    delta=1/1000
    for it in range(1001):
      F=beta*np.matmul(C,[I_V,I_N])
      T=[S_V*F[0],S_N*F[1]]
      I_V+=delta*(T[0]-gamma*I_V)
      I_N+=delta*(T[1]-gamma*I_N)
      S_V+=-delta*(T[0]+alpha)
      S_N+=-delta*T[1]
      z=log(I_V/I_V0)/(delta*(it+1))
      print(z,(z-L_V1)/(delta*(it+2))*2)
    print()
    
