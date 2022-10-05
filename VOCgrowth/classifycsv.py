import sys,argparse,csv

parser=argparse.ArgumentParser()
parser.add_argument('-b', '--maxbad',      type=float, default=0.05, help="Maximum proportion of Ns allowed")
parser.add_argument('-c', '--compare',     action="store_true",      help="Compare with original lineage")
parser.add_argument('-f', '--mindate',     default="2019-01-01",     help="Min sample date of sequence")
parser.add_argument('-g', '--gisaid',      action="store_true",      help="Use GISAID data instead of COG-UK data")
parser.add_argument('-t', '--maxdate',     default="9999-12-31",     help="Max sample date of sequence")
parser.add_argument('-s', '--sorted',      action="store_true",      help="Guarantees the input is sorted in reverse date order")
args=parser.parse_args()

from classify import classify, expandlin, contractlin
  
def mutationlist(mutations):
  if args.gisaid: return mutations[1:-1].split(',')
  else: return mutations.split('|')

def badness(Ncontent):
  if args.gisaid:
    if Ncontent!="": return float(Ncontent)
    else: return 0
  else:
    t=0
    for x in Ncontent.split('|'):
      y=x.split('-')
      if len(y)==1: t+=1
      else: t+=int(y[1])-int(y[0])+1
    return t/29903

if args.gisaid:
  sep='\t'
  keys=["Collection date","Pango lineage","AA Substitutions","N-Content"]
else:
  sep=','
  keys=["sample_date","lineage","mutations","ambiguities"]
t0=t1=0
reader=csv.reader(sys.stdin,delimiter=sep)
writer=csv.writer(sys.stdout,delimiter=sep)
headings=next(reader)
writer.writerow(headings)
dec={h:i for (i,h) in enumerate(headings)}
cols=[dec[h] for h in keys]

comp={}
while 1:
  row=next(reader,None)
  if row==None: break
  (date,lineage,mutations,Ncontent)=[row[i] for i in cols]
  if len(date)<10: continue
  if date>args.maxdate: continue
  if date<args.mindate:
    if args.sorted: break
    continue
  bad=badness(Ncontent)
  t0+=1
  if bad>args.maxbad: t1+=1;continue
  newlin=classify(mutationlist(mutations),lineage)
  row[cols[1]]=newlin
  if args.compare:
    if lineage not in comp: comp[lineage]=[0]*len(lineages)
    comp[lineage][lineages.index(newlin)]+=1
  writer.writerow(row)

print("Discarded",t1,"from",t0,"(%.1f%%) due to bad coverage"%(t1/t0*100),file=sys.stderr)

if args.compare:
  # Wildcard ending is replaced with '.'. It's a match if it's equal to a prefix of (database lineage)+'.'
  # Note that BA.5* will match BA.5.1 and BA.5 but not BA.53
  #           BA.5.* will match BA.5.1 but not BA.5 or BA.53
  targlinsexact=[]
  targlinsprefix=[]
  for lin in lineages:
    if lin[-2:]=='.*': exact="-";prefix=lin[:-1]
    elif lin[-1]=='*': exact=lin[:-1];prefix=lin[:-1]+'.'
    else: exact=lin;prefix="-"
    targlinsexact.append(expandlin(exact))
    targlinsprefix.append(expandlin(prefix))
  
  l=[(sum(comp[oldlin]),oldlin) for oldlin in comp]
  l.sort(reverse=True)
  # Unassigned
  match=nonmatch=0
  wid=[max(len(lin)+1,7) for lin in lineages]
  mwid=max(max(wid),8)
  dbwid=max(len(oldlin) for oldlin in comp)
  print(file=sys.stderr)
  print("%*s  %*s  "%(dbwid,"DB_lin",mwid,"DB_group")," ".join("%*s"%(w,lin) for (w,lin) in zip(wid,lineages)),"     Match  Nonmatch",file=sys.stderr)
  rest=0
  for (tot,oldlin) in l:
    if oldlin=="Unassigned":
      if tot>=10: print("%*s"%(dbwid,oldlin)," %*s"%(mwid,"[any]")," "," ".join("%*s"%(w,n) for (w,n) in zip(wid,comp[oldlin]))," %9s %9s"%("-","-"),file=sys.stderr)
      else: rest+=1
    else:
      dblin=expandlin(oldlin)
      ind=len(lineages)-1
      for i in range(len(lineages)):
        exact=targlinsexact[i]
        prefix=targlinsprefix[i]
        if dblin==exact:
          if prefix=='-': ind=i;break# Exact match with non-wildcard takes precedence over anything later
        if (dblin+'.')[:len(prefix)]==prefix: ind=i
      good=comp[oldlin][ind]
      bad=sum(comp[oldlin])-comp[oldlin][ind]
      if tot>=10: print("%*s"%(dbwid,oldlin)," %*s"%(mwid,lineages[ind])," "," ".join("%*s"%(w,n) for (w,n) in zip(wid,comp[oldlin]))," %9d %9d"%(good,bad),file=sys.stderr)
      else: rest+=1
      match+=good
      nonmatch+=bad
  if rest>0:
    print("+ %d other%s"%(rest,"s" if rest!=1 else ""),file=sys.stderr)
  print("Total match: %d (%.2f%%)"%(match,match/(match+nonmatch)*100),file=sys.stderr)
  print("Total nonmatch: %d (%.2f%%)"%(nonmatch,nonmatch/(match+nonmatch)*100),file=sys.stderr)
