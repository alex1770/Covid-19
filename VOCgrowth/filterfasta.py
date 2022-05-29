# Picks out genomes whose ID is in idlist
# The ID is defined by the bit before the '|', if such exists.
# IDs in idlist can include the initial '>' or not.
# 
# python filterfasta.py idlist < input.fasta > output.fasta

import sys

with open(sys.argv[1]) as fp:
  okids=set(x.lstrip('>') for x in fp.read().strip().split('\n'))

# mode 0: skipping lines until the next '>'
#      1: outputting lines until the next '>'
mode=0
for x in sys.stdin:
  if x[0]=='>':
    f=x.find('|')
    if f>=0: id=x[1:f]
    else: id=x.rstrip()[1:]
    if id in okids:
      print(x,end='')
      mode=1
    else:
      mode=0
  else:
    if mode: print(x,end='')
