from subprocess import Popen,PIPE

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): p.write((' '.join(map(str,s))+'\n').encode('utf-8'))

fn="hitout"
with open(fn) as fp:
  ll=[l.split() for l in fp.read().strip().split('\n')]
R0=ll[0][3]
distdesc=ll[1][1:]

for mode in [0,1]:
  desc=["susceptibility","connectivity"][mode]
  outfn="HITgraph_"+desc+".png"
  p=Popen("gnuplot",shell=True,stdin=PIPE).stdin
  write('set terminal pngcairo font "sans,13" size 1920,1280')
  write('set bmargin 5;set lmargin 15;set rmargin 15;set tmargin 5')
  write('set output "%s"'%outfn)
  write('set xtics nomirror')
  write('set y2tics mirror')
  title='HIT vs CV in %s mode at R_0 = %s'%(desc,R0)
  write('set title "%s"'%title)
  write('set xlabel "Coefficient of Variation"')
  write('set ylabel "HIT"')
  #write('set grid ytics lc rgb "#dddddd" lt 1')
  s='plot '
  #write('plot "-" using 1:xtic(1) with lines')
  for dist in range(4):
    if s!='plot ': s+=', '
    s+='"-" using 1:2 with lines lw 2 title "%s"'%distdesc[1+mode+dist*2]
  write(s)
  for dist in range(4):
    for l in ll:
      if l[0]!='#': write(l[0],l[1+mode+dist*2])
    write('e')
  p.close()
  print("Written graph to %s"%outfn)
