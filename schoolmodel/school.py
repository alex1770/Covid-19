from math import log,exp,sqrt
import numpy as np

#                                   Child    Adult
# Basic contact matrix C = Child [  C_cc     C_ca   ]
#                          Adult [  C_ac     C_aa   ]
#
#
# N_i = number of people in group i (i=c or a)
# C_{ij} = (number of people from group j encountered by an individual in group i)
# C_ca/C_ac = N_a/N_c
#
# I_c1 = children infected by children
# I_c2 = children infected by adults
# I_c1, I_c2, S_c, I_a, S_a represent absolute numbers of their respective populations
#
# dI_c1/dt = S_c*C_cc*(beta_cc1.I_c1+beta_cc2.I_c2)/N_c
# dI_c2/dt = S_c*C_ca*(beta_ca.I_a)/N_a
# dI_a/dt =  S_a*(C_ac*(beta_ac1.I_c1+beta_ac2.I_c2)/N_c+C_aa*beta_aa.I_a/N_a)
#
# beta_ij = P(person from j is not isolating|infected)*P(encounter would cause an infection j->i | person from j is infected and not isolating)
#         = P(person from j is not isolating|infected)*Transmissibility(j)*Susceptibility(i)
# (i=c,a; j=c1,c2,a)
#
# Leads to derived contact matrix
#
#                     c1              c2             a
#          D =  c1 [  C_cc.beta_cc1   C_cc.beta_cc2                ]
#               c2 [                                 C_ca.beta_ca  ]
#                a [  C_ac.beta_ac1   C_ac.beta_ac2  C_aa.beta_aa  ]
#
# so that
#
#      [ F_c1 ] = D . [ I_c1/N_c ]
#      [ F_c2 ]       [ I_c2/N_c ]
#      [ F_a  ]       [ I_a/N_a  ]
#
# and
#
# d/dt [ I_c1 ]  =  [ S_c.F_c1 ] - gamma.[ I_c1 ]
#      [ I_c2 ]     [ S_c.F_c2 ]         [ I_c2 ]
#      [ I_a  ]     [ S_a.F_a  ]         [ I_a  ]
#
# d/dt [ S_c ]  =  - [ S_c.(F_c1+F_c2) ]
#      [ S_a ]       [ S_a.F_a         ]

# Order of indices is c,a or c1,c2,a

######################################################################################
# PARAMETERS

# https://explore-education-statistics.service.gov.uk/find-statistics/school-pupils-and-their-characteristics
N = [ 8.9e6, 57.8e6 ]

# Guesstimates (not to be relied on)

suscep = np.array([ 0.25, 0.25 ])
transm = np.array([ 0.15, 0.15, 0.25 ])
nonisolate = [ 0.5, 0.2, 0.5]# Probability of not isolating given infected

C = np.array([[ 8,  3],
              [ -1, 3]], dtype=float)
C[1,0]=C[0,1]*N[0]/N[1]

beta=np.zeros((2,3))
for i in range(2):
  for j in range(3):
    beta[i,j]=transm[j]*suscep[i]*nonisolate[j]

# ONS infection survey, 18 Sep shows ~0.3% prevalence
I = np.array([ 0.0015*N[0], 0.0015*N[0], 0.003*N[1]])

S = np.array([ 0.85*N[0], 0.85*N[1]])

gamma = 0.1 # recovery rate

#######################################################################################

D = np.array([[C[0,0]*beta[0,0], C[0,0]*beta[0,1], 0],
              [0,                0,                C[0,1]*beta[0,2]],
              [C[1,0]*beta[1,0], C[1,0]*beta[1,1], C[1,1]*beta[1,2]]])

NN=np.array([N[0], N[0], N[1]])

print('"Child1" means Children infected by children')
print('"Child2" means Children infected by adults')
print()
print("Susceptibilities (child, adult):",suscep)
print("Transmissibilities (child1, child2, adult):",transm)
print("Non-isolation probabilities (child1, child2, adult):",nonisolate)
print()
print("Simple contact matrix:")
print(C)
print()
print("Derived contact matrix:")
print(D)
print()

print("I_c1 = Children infected by children, as a proportion of all children")
print("I_c1 = Children infected by adults, as a proportion of all children")
print("I_a  = Infected adults, as a proportion of all adults")
print()
print(" Day   I_c1 %%       I_c2 %%       I_a %%")

subdiv=1000
delta=1/subdiv
days=60
for step in range(days*subdiv+1):
  if step%subdiv==0: print("%4d"%(step//subdiv),'   '.join("%10.6f"%(x*100) for x in I/NN))
  F=np.matmul(D,I/NN)
  T=[S[0], S[0], S[1]]*F
  I+=delta*(T-gamma*I)
  S-=delta*np.array([T[0]+T[1], T[2]])
