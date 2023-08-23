s=""
with open("nucleotide_counts") as fp:
  for x in fp:
    y=x.split()
    l=[]
    for z in y[1:]:
      c=z[0]
      n=int(z[1:])
      l.append((n,c))
    (n,c)=max(l)
    s+=c
#print(s)

with open("refgenome") as fp:
  ref=fp.read()

n=len(ref)
assert len(s)==n
for i in range(n):
  if s[i]!=ref[i]: print(i,ref[i],s[i])
  
