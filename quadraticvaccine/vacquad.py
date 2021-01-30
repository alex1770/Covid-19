from math import log,exp
import numpy as np

# Two population (V=to-be-vaccinated, N=not-to-be-vaccinated) SIR with vaccine term (alpha)
# Assume a simple interaction matrix beta*[[p,1-p],[1-p,p]]:
#
# (d/dt) [ I_V ] = ( beta.[ pS_V   (1-p)S_V ] - gamma.[ 1  0 ] ).[ I_V ]
#        [ I_N ]          [ (1-p)S_N   pS_N ]         [ 0  1 ]   [ I_N ]
#  
# (d/dt) [ S_V ] = - beta.[ pS_V   (1-p)S_V ].[ I_V ] - alpha.[ 1 ]
#        [ S_N ]          [ (1-p)S_N   pS_N ] [ I_N ]         [ 0 ]
#
# Aim to work out I_V, I_N to second order in t in order to get expressions of the form
# I_V(t) = I_V0*exp(L_V1*t+(1/2)L_V2*t^2)
# I_N(t) = I_N0*exp(L_N1*t+(1/2)L_N2*t^2)


I_N0=0.02
I_V0=0.01
S_N0=0.8
S_V0=0.9
beta=0.09
gamma=0.09
alpha=0.006*67/3

print("#    p           L_V1          L_V2          L_N1          L_N2")

for p in [0.1*x for x in range(11)]:
  
  # First work out S_V, S_N to first order in t
  # S_V = S_V0 + S_V1.t + ...
  # S_N = S_N0 + S_N1.t + ...
  A=np.array([[p*S_V0, (1-p)*S_V0],
             [(1-p)*S_N0, p*S_N0]])
  [S_V1, S_N1] = -beta*np.matmul(A,[I_V0,I_N0]) - np.array([alpha,0])
  
  # (d/dt) [I_V] = ( beta.[ pS_V   (1-p)S_V ] - gamma.[ 1  0 ] ).[I_V] =  (B + Ct).[I_V] + O(t^2)
  #        [I_N]          [ (1-p)S_N   pS_N ]         [ 0  1 ]   [I_N]             [I_N]
  #
  B=beta*A - gamma*np.array([[1,0],[0,1]])
  C=beta*np.array([[p*S_V1, (1-p)*S_V1],
                  [(1-p)*S_N1, p*S_N1]])

  # Solving above DE gives:
  # (I_V,I_N)^T = (I + Bt + (1/2)(B^2+C)t^2 + ...) (I_V0, I_N0)^T
  # 
  # I_V = I_V0 + I_V1.t + (1/2)I_V2.t^2 + ...
  # I_N = I_N0 + I_N1.t + (1/2)I_N2.t^2 + ...
  #
  [I_V1, I_N1] = np.matmul(B,[I_V0,I_N0])
  [I_V2, I_N2] = np.matmul(np.matmul(B,B)+C,[I_V0,I_N0])

  # log(I_V) = log(I_V0) + I_V1/I_V0.t + (1/2)(I_V2/I_V0-(I_V1/I_V0)^2).t^2 + ...
  # log(I_N) = log(I_N0) + I_N1/I_N0.t + (1/2)(I_N2/I_N0-(I_N1/I_N0)^2).t^2 + ...
  #
  L_V0=log(I_V0);L_V1=I_V1/I_V0;L_V2=I_V2/I_V0-(I_V1/I_V0)**2
  L_N0=log(I_N0);L_N1=I_N1/I_N0;L_N2=I_N2/I_N0-(I_N1/I_N0)**2
  print("%6.3f   %12.8f  %12.8f  %12.8f  %12.8f"%(p,L_V1,L_V2,L_N1,L_N2))
  if 0:
    print("B:")
    print(B)
    print()
    print("C:")
    print(C)
    print()
    print("B^2-C:")
    print(np.matmul(B,B)-C)
    print()
    print("--------------------")
    print()
  
