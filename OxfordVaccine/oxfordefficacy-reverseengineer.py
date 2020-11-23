e0=0.90
e1=0.62
e01=0.70
n=131

r0=1-e0
r1=1-e1
r01=1-e01

p0=n*(r1-r01)/((r1-r0)*(1+r01))
p1=n*(r01-r0)/((r1-r0)*(1+r01))
v0=r0*p0
v1=r1*p1

p0=int(p0+.5)
p1=int(p1+.5)
v0=int(v0+.5)
v1=int(v1+.5)
print(p0,v0,p1,v1)
print(p0+p1+v0+v1)
print(1-v0/p0)
print(1-v1/p1)
print(1-(v0+v1)/(p0+p1))

    
for p0 in range(25,35):
  for v0 in range(10):
    for p1 in range(60,80):
      v1=n-(p0+v0+p1)
      if abs(v0/p0-r0)<0.00501 and abs(v1/p1-r1)<0.00501 and abs((v0+v1)/(p0+p1)-r01)<0.00501:
        print(p0,v0,p1,v1)
