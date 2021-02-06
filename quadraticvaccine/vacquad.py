from math import log,exp,sqrt
import numpy as np

# Two population (V=to-be-vaccinated, N=not-to-be-vaccinated) SIR with vaccine term (alpha)
# Assume a simple contact matrix C = k.[p  1-p].[v     ]
#                                      [1-p  p] [   1-v]
#
# where v is the proportion of people in the vaccinated group, and k is a scale factor, 1/(p/2+sqrt(p^2/4-(2p-1)v(1-v))), chosen so that
# k.[p  1-p].[v    0] has its largest eigenvalue equal to 1, so that R_0=beta/gamma as usual.
#   [1-p  p] [0  1-v]
#
# I_V, S_V, I_N, S_N represent proportions of their respective populations (so between 0 and 1)
#
# Force term:
#        [ F_V ] = beta.k.[ p   1-p ].[v     ].[ I_V ] = beta.C.[ I_V ]
#        [ F_N ]          [ 1-p   p ] [   1-v] [ I_N ]          [ I_N ]
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
alpha=0.006/v

print("#    p           L_V1          L_V2          L_N1          L_N2")

#for p in [0.02*x for x in range(25,51)]:
for p in [0.9]:

  C=np.array([[p*v,(1-p)*(1-v)],[(1-p)*v,p*(1-v)]])
  C*=1/(p/2+sqrt(p**2/4-(2*p-1)*v*(1-v)))
  
  # First work out S_V, S_N to first order in t
  # S_V = S_V0 + S_V1.t + ...
  # S_N = S_N0 + S_N1.t + ...
  F0=beta*np.matmul(C,[I_V0,I_N0])
  [S_V1, S_N1] = -([S_V0,S_N0]*F0) - np.array([alpha,0])
  
  # (d/dt) [I_V] = ( beta.[ S_V     ].C - gamma.[ 1    ] ).[I_V] = (say) (A + Bt + O(t^2)).[I_V]
  #        [I_N]          [     S_N ]           [    1 ]   [I_N]                           [I_N]
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

  if 1:
    (I_V,S_V,I_N,S_N)=(I_V0,S_V0,I_N0,S_N0)
    delta=1/1000
    for it in range(1001):
      F=beta*np.matmul(C,[I_V,I_N])
      T=[S_V*F[0],S_N*F[1]]
      I_V+=delta*(T[0]-gamma*I_V)
      I_N+=delta*(T[1]-gamma*I_N)
      S_V+=-delta*(T[0]+alpha)
      S_N+=-delta*T[1]
      print(log(I_V/I_V0)/(delta*(it+1)))
      if it==50: x=log(I_V/I_V0)
      if it==100: y=log(I_V/I_V0)
    print("Check:",10*(x*4-y),400*(y-2*x))
    print()
    
