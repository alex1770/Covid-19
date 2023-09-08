# Applies a base score, such as derived from intersectional_find_rare_mutations.py,
# to each aligned fasta entry.

import sys
from math import log

# These score BA.2.86-ness
# Base location counting from 1
scorestext="""
15756 A    87
22556 G   134
  897 A   138
26610 G   160
22896 A   260
 7842 G   316
28958 A   320
22353 A   363
22916 T  1027
21711 T  1155
21640 A  1221
23005 A  1251
26529 C  1354
21639 A  1456
21941 T  1807
"""

numg=15e6# Approx number of genomes in GISAID at the time nucleotide_counts was made
scores=[]
for x in scorestext.strip().split('\n'):
  y=x.split()
  scores.append([int(y[0])-1,y[1],log(numg/int(y[2]))])

for x in sys.stdin:
  assert x[0]=='>'
  f=x.find('|')
  if f>=0: id=x[1:f]
  else: id=x.rstrip()[1:]
  seq=next(sys.stdin).strip()
  sc=0
  for (p,b,v) in scores:
    if seq[p]==b: sc+=v
  nn=len([x for x in seq if x=='N' or x=='-'])
  print("%40s   %7.2f   %5d"%(id,sc,nn))
