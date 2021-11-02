from stuff import *

genes={
  'ORF1a': [266, 13467],
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
  'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
  'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}
  
acgt={}
with open('aligned.examples.fasta') as fp:
  for id in fp:
    id=id[1:].strip()
    seq=next(fp)
    acgt[id]=[seq,None]

with open('refgenome') as fp:
  refgenome=fp.read()

me=loadcsv("mutationexamples.csv")
for (id,ml) in zip(me['sequence_name'],me['mutations']):
  acgt[id][1]=ml.split('|')

for id in acgt:
  (seq,ml)=acgt[id]
  for mut in ml:
    (gene,m)=mut.split(':')
    x=m[0]
    p1=int(m[1:-1])
    y=m[-1]
    if gene=='synSNP':
      p=p1-1
      continue
    else:
      if gene=='orf1ab':
        if p1<=4401: gene='ORF1a'
        else: gene='ORF1b';p1-=4401
      p0=genes[gene]
      p=p0[0]+(p1-1)*3-1
      continue
    from_acgt=codontable.get(refgenome[p:p+3],'-');from_mut=x
    to_acgt=codontable.get(seq[p:p+3],'-');to_mut=y
    print(from_acgt,from_mut,to_acgt,to_mut)
    assert from_acgt==from_mut and to_acgt==to_mut
