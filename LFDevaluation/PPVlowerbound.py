import os,json

lfddir='../apidata_lfd'

# Looks very much like the erroneous 2021-03-26 file (published 2021-03-27) contains the true full LFDcases
#captureddates=sorted(os.listdir(lfddir))
captureddates=sorted(x for x in os.listdir(lfddir) if x[:2]=='20')

dd={}
for fn in captureddates:
  with open(os.path.join(lfddir,fn),'r') as fp:
    dd[fn]={y['date']: y for y in json.load(fp)}

# dd[publishdate][historicaldate]={historical data}
    
PPVLB={date:0 for date in captureddates}
n=len(captureddates)
for i0 in range(n-1):
  date0=captureddates[i0]
  for i1 in range(i0+1,n):
    date1=captureddates[i1]
    for date2 in captureddates:
      if date2 in dd[date0] and date2 in dd[date1]:
        ppv=1-dd[date1][date2]['newCasesLFDOnlyBySpecimenDate']/dd[date0][date2]['newCasesLFDOnlyBySpecimenDate']
        if ppv>PPVLB[date2]: PPVLB[date2]=ppv

d=dd[captureddates[-1]]

for x in sorted(list(PPVLB)):
  a=d[x]['newCasesLFDOnlyBySpecimenDate']
  b=d[x]['newCasesLFDConfirmedPCRBySpecimenDate']
  print(x,"PPV >= %5.1f%% %5.1f%%"%(PPVLB[x]*100,b/(a+b)*100))
