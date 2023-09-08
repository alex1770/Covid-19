# Convert COG metadata into possible (aligned) fasta that would produce these mutations and ambiguities
# It will choose some plausible triple of bases for each AA name, which will often not be the actual triple used.

from stuff import *
import sys,argparse

parser=argparse.ArgumentParser()
parser.add_argument('-f', '--mindate',       default="2019-01-01",  help="Minimum date used from input")
parser.add_argument('-s', '--sorted',     action="store_true",      help="Whether input is sorted into reverse date order")
args=parser.parse_args()


# Counting from 1, inclusive-exclusive notation
genes={
  'ORF1a': [266, 13469],
  'ORF1b': [13468, 21556],
  'S': [21563, 25385],
  'ORF3a': [25393, 26221],
  'E': [26245, 26473],
  'M': [26523, 27192],
  'ORF6': [27202, 27388],
  'ORF7a': [27394, 27760],
  'ORF7b': [27756, 27888],
  'ORF8': [27894, 28260],
  'N': [28274, 29534],
  'ORF10': [29558, 29675]
}
# ORF1ab codon position p>4401 <-> ORF1b codon position p-4401

codontable = {
  'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
  'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
  'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

with open('refgenome') as fp:
  refgenome=fp.read()

stats={}
for geneloc in genes.values():
  for p in range(geneloc[0],geneloc[1],3):
    c=refgenome[p:p+3]
    x=codontable.get(c,'?')
    if x not in stats: stats[x]={}
    stats[x][c]=stats[x].get(c,0)+1


# Find most likely triple of bases for each AA
# AA2bases[ new AA name, refgenome triple ] = best guess for new triple
AA2bases={}
w0,w1,w2=0.3,0.01,0.3# Weights that work reasonably well
for x in stats:
  y=stats[x]
  best=(-1e30,)
  s=sum(stats[x].values())
  for old3 in codontable:
    for new3 in y:
      # Mutation from old3 -> new3
      score=y[new3]/s-w0*(old3[0]!=new3[0])-w1*(old3[1]!=new3[1])-w2*(old3[2]!=new3[2])
      if score>best[0]: best=(score,new3)
    AA2bases[x,old3]=best[1]

for (name,date,ml,am) in csvrows_it(sys.stdin,["sequence_name","sample_date","mutations","ambiguities"]):
  if date<args.mindate:
    if args.sorted: break
    continue
  
  ml=ml.split('|')
  am=am.split('|')

  output=list(refgenome)
  
  for mut in ml:
    (gene,m)=mut.split(':')
    x=m[0]
    p1=int(m[1:-1])
    y=m[-1]
    if gene=='synSNP':
      p=p1-1
      output[p]=y
    else:
      if gene=='orf1ab':
        if p1<=4401: gene='ORF1a'
        else: gene='ORF1b';p1-=4401
      p0=genes[gene]
      p=p0[0]+(p1-1)*3-1
      output[p:p+3]=list(AA2bases[y,refgenome[p:p+3]])

  for ar in am:
    x=[int(p)-1 for p in ar.split('-')]
    if len(x)==1: output[x[0]]='N'
    else: output[x[0]:x[1]+1]=['N']*(x[1]+1-x[0])
    
  print(">hCoV-19/"+name)
  print("".join(output))
