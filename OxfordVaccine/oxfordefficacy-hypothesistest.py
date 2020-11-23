
(p0,v0,p1,v1)=(30,3,71,27)

r=(v0+v1)/(p0+p1)
print("H0: Efficacy = %.1f%% regardless of which of the two dose regimens are used"%((1-r)*100))

from scipy.stats import poisson, binom

# Assume p0 was really from a Poisson(p0), and assume v0 was from a B(Poisson(p0), r)
# Assume p1 was really from a Poisson(p1), and assume v1 was from a B(Poisson(p1), r)
N=1000000
P0=poisson.rvs(p0,size=N)
V0=binom.rvs(poisson.rvs(p0,size=N),r)
P1=poisson.rvs(p1,size=N)
V1=binom.rvs(poisson.rvs(p1,size=N),r)

delta=v1/p1-v0/p0
n=(V1/P1-V0/P0>=delta).sum()
print("Chance of difference of efficacies being as observed or greater under H0 is %.2f%%"%(n/N*100))

delta=v1/p1/(v0/p0)
n=(V1/P1>=delta*V0/P0).sum()
print("Chance of ratio of efficacies being as observed or greater under H0 is %.2f%%"%(n/N*100))

