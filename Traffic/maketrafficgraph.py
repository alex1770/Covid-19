from subprocess import Popen,PIPE

# Use this to cater for earlier versions of Python whose Popen()s don't have the 'encoding' keyword
def write(*s): inp.write((' '.join(map(str,s))+'\n').encode('utf-8'))

outfn="trafficgraph.png"
p=Popen("gnuplot",shell=True,stdin=PIPE);inp=p.stdin
write('set terminal pngcairo font "sans,10" size 2400,640')
write('set bmargin 5;set lmargin 12;set rmargin 12;set tmargin 4')
write('set output "%s"'%outfn)
write('set tics nomirror')
title='Congestion levels in London as derived from Google maps images; vertical lines are at the start of Mondays (midnight)'
write('set title "%s"'%title)
#write('set xlabel "blah1"')
#write('set ylabel "blah2"')
write('set key left')
write('set timefmt "%Y-%m-%dT%H:%M:%S"')
write('set xdata time')
write('set format x "%Y-%m-%d"')
write('set xrange ["2020-03-20":]')
write('set tics scale 2,0.5')
write('set grid xtics noytics lc rgb "#dddddd" lt 1')
write('set xtics "2020-03-16", 604800')
write('plot "trafficlevels" u 2:($5)+3*($6)+10*($7) w lines')
write('quit')
inp.close()
p.wait()
print("Written graph to %s"%outfn)
