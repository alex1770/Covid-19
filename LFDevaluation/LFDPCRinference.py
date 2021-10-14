import os,json,sys
from scipy.linalg import block_diag
from math import log
import numpy as np
from stuff import *

np.set_printoptions(precision=3,suppress=True)
np.set_printoptions(edgeitems=30, linewidth=100000)

# Try to reverse-engineer PCR-retest-negatives-by-specimen-date using the assumption that the LFD+ves and PCR-retest+ves follow newLFDTests
# in a uniform way, and cross-checking with PCR-retest-negatives-by-publication-date

# https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaName=England&metric=newCasesLFDConfirmedPCRBySpecimenDate&metric=newCasesLFDOnlyBySpecimenDate&metric=newCasesPCROnlyBySpecimenDate&metric=newLFDTests&metric=newCasesBySpecimenDate&format=csv&release=2021-04-18
# https://docs.google.com/spreadsheets/d/1yItIPgFTItAM29BZ6D8unah53EkOWMkiCW-T5TcGePI/edit#gid=0

lfddir='../apidata_lfd'

#mindate='2021-05-18'
mindate='2021-08-16'
#mindate='2021-09-28'

dates=sorted(x for x in os.listdir(lfddir) if x[:2]=='20' and x>=mindate)

dd={}
for fn in dates:
  with open(os.path.join(lfddir,fn),'r') as fp:
    dd[fn]={y['date']: y for y in json.load(fp) if y['date']>=mindate}

# dd[adjustedpublishdate][historical specimen date]={historical data}

now=dates[-1]

# nhist = number of days of LFD history to use to predict retest results
if len(sys.argv)>1: nhist=int(sys.argv[1])
else: nhist=2

# weight = how strongly to force results-by-specimen-date to agree with results-by-publication-date
if len(sys.argv)>2: weight=float(sys.argv[2])
else: weight=0

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  for t in range(3):
    response = requests.get(url+req, timeout=30)
    if response.ok: break
  else: raise RuntimeError('Request failed: '+response.text)
  return response.json()['body'][::-1]

d=datetime.datetime.now(pytz.timezone("Europe/London"))
today=datetoday(d.strftime('%Y-%m-%d'))
if d.hour+d.minute/60<16+15/60: today-=1# Dashboard/api updates at 4pm UK time; allow 15 minutes for things to settle down
minday=datetoday(mindate)

cachedir_casesbypubdate='../apidata_casesbypubdate'
os.makedirs(cachedir_casesbypubdate,exist_ok=True)
todate=daytodate(today)
fn=os.path.join(cachedir_casesbypubdate,todate)
if os.path.isfile(fn):
  with open(fn,'r') as fp: cases=json.load(fp)
else:
  cases=get_data('areaType=nation&areaCode=E92000001&metric=newCasesByPublishDate&metric=cumCasesByPublishDate&format=json')
  with open(fn,'w') as fp: json.dump(cases,fp,indent=2)
  print("Retrieved api data at publication date",todate)

cumcasespub=np.zeros(today-minday+1,dtype=int)
newcasespub=np.zeros(today-minday+1,dtype=int)
for x in cases:
  i=datetoday(x['date'])-minday
  if i>=0:
    cumcasespub[i]=x['cumCasesByPublishDate']
    newcasespub[i]=x['newCasesByPublishDate']

removed=newcasespub[1:]-(cumcasespub[1:]-cumcasespub[:-1])

l=[]
for day in range(minday-1,today):
  fn=os.path.join(lfddir,daytodate(day))
  with open(fn,'r') as fp:
    l.append(sum(x['newCasesLFDConfirmedPCRBySpecimenDate'] for x in json.load(fp)))
lfdconfirmed=np.array(l)
lfdnewconfbypub=lfdconfirmed[1:]-lfdconfirmed[:-1]
for i in range(today-minday-1):
  if lfdnewconfbypub[i]<=0: t=lfdnewconfbypub[i]+lfdnewconfbypub[i+1];lfdnewconfbypub[i]=t//2;lfdnewconfbypub[i+1]=(t+1)//2

if 0:
  maxdays=10
  for specdate in dates[:-3]:
    nn=np.array([dd[repdate][specdate]['newLFDTests'] for repdate in dates if repdate>=specdate])
    oo=np.array([dd[repdate][specdate]['newCasesLFDOnlyBySpecimenDate'] for repdate in dates if repdate>=specdate])
    cc=np.array([dd[repdate][specdate]['newCasesLFDConfirmedPCRBySpecimenDate'] for repdate in dates if repdate>=specdate])
    nn=nn[:maxdays]
    oo=oo[:maxdays]
    cc=cc[:maxdays]
    # Try to find (e.g., nhist=3) best a,b,c,lam s.t. a*nn+b*(0,nn[:-1])+c*(0,0,nn[:-2]) is close to oo+lam*cc
    # lam should be 1/PPV
    n=len(nn)
    specday=datetoday(specdate)
    badrows=[i for i in range(n) if removed[specday+i-minday]==0 or specday+i==datetoday('2021-10-01')]
    nn=np.delete(nn,badrows)
    oo=np.delete(oo,badrows)
    cc=np.delete(cc,badrows)
    n=len(nn)
    nnl=[]
    for i in range(nhist):
      nnl.append(np.concatenate([[0]*i,nn[:len(nn)-i]]))
    nnl.append(-cc)
    nnn=np.transpose(np.array(nnl,dtype=float))
    dnnn=nnn-np.concatenate([[[0]*(nhist+1)],nnn])[:n]
    doo=oo-np.concatenate([[0],oo])[:n]
    condition=100
    dnnn[:,:nhist]=dnnn[:,:nhist]/condition
    rr=np.linalg.lstsq(dnnn,doo)
    r=rr[0];r[:nhist]=r[:nhist]/condition
    print(specdate,"%5.3f %6.1f"%(1/r[nhist],cc[-1]*(r[nhist]-1)),end='')
    ccn=cc/nn;oon=oo/nn
    for i in range(min(len(nn)-1,5)):
      ppv=(ccn[-1]-ccn[i])/(oon[i]-oon[-1])
      print(" %5.3f"%ppv,end='')
    print(" ",rr[1],r[:nhist]*1000,r[1:nhist]/r[0]*100)
    
if 1:
  stabledays=10# Number of days after which you don't expect the specimen day metrics to change
  tnnn=too=None
  lastday=today-nhist
  nspec=lastday-minday
  baddays=set(i for i in range(today-minday) if removed[i]==0)
  baddays.add(datetoday('2021-10-01')-minday)
  eqq=np.zeros([nspec-stabledays,nspec*(nhist+1)])
  removedadj=removed.copy()
  for i in baddays: removedadj[i+1]+=removed[i];removedadj[i]=0
  eqo=removedadj[stabledays:nspec].copy()
  for specday in range(minday,lastday):
    # i^th element of arrays corresponds to publiction day specday+i
    specdate=daytodate(specday)
    s=specday-minday# specimen day (as an offset from minday)
    nn=np.array([dd[repdate][specdate]['newLFDTests'] for repdate in dates if repdate>=specdate])
    oo=np.array([dd[repdate][specdate]['newCasesLFDOnlyBySpecimenDate'] for repdate in dates if repdate>=specdate])
    cc=np.array([dd[repdate][specdate]['newCasesLFDConfirmedPCRBySpecimenDate'] for repdate in dates if repdate>=specdate])
    drr=removed[s:][:stabledays]
    nn=nn[:stabledays]
    oo=oo[:stabledays]
    cc=cc[:stabledays]
    pp=np.arange(s,s+stabledays)# list of publication days (as an offset from minday)
    badrows=[i for i in range(len(drr)) if s+i in baddays]
    # Remove bad reporting days
    nn=np.delete(nn,badrows)
    oo=np.delete(oo,badrows)
    cc=np.delete(cc,badrows)
    pp=np.delete(pp,badrows)
    drr=np.delete(drr,badrows)
    # Make delta by publication date
    n=len(nn)
    dnn=nn-np.concatenate([[0],nn])[:n]
    doo=oo-np.concatenate([[0],oo])[:n]
    dcc=cc-np.concatenate([[0],cc])[:n]
    # Try to find (e.g., nhist=3) best a,b,c,lam s.t. a*nn+b*(0,nn[:-1])+c*(0,0,nn[:-2]) is close to oo+lam*cc
    # lam should be 1/PPV
    dnnl=[]
    for i in range(nhist):
      dnnl.append(np.concatenate([[0]*i,dnn[:n-i]]))
    dnnl.append(-dcc)
    dnnn=np.transpose(np.array(dnnl,dtype=float))
    condition=100
    dnnn[:,:nhist]=dnnn[:,:nhist]/condition
    #rr=np.linalg.lstsq(dnnn,doo)
    if tnnn is None:
      tnnn=dnnn;too=doo
    else:
      tnnn=block_diag(tnnn,dnnn)
      too=np.concatenate([too,doo])
    for i in range(n):
      s=specday-minday# s=spec day counted from minday
      p=pp[i]         # p=pub day counted from minday
      if p>=stabledays and p<nspec:
        eqq[p-stabledays,s*(nhist+1)+nhist]+=dcc[i]
        eqo[p-stabledays]+=dcc[i]
  tnnn=np.concatenate([tnnn,eqq*weight])
  too=np.concatenate([too,eqo*weight])
  trr=np.linalg.lstsq(tnnn,too)
  tr=trr[0]
  #z=eqq@tr
  #for x in zip(dates[stabledays:-3],eqo,z): print(x)
  removedadj2=removed.copy()
  for i in baddays: a=(removed[i]+removed[i+1])/2;removedadj2[i]=a;removedadj2[i+1]=a
  ccpub=np.array([sum(dd[repdate][specdate]['newCasesLFDConfirmedPCRBySpecimenDate'] for specdate in dates if repdate>=specdate) for repdate in dates])
  dccpub=ccpub-np.concatenate([[0],ccpub])[:len(ccpub)]
  for day in range(minday,lastday):
    specdate=daytodate(day)
    s=day-minday
    r=tr[s*(nhist+1):(s+1)*(nhist+1)]
    r[:nhist]=r[:nhist]/condition
    nn=np.array([dd[repdate][specdate]['newLFDTests'] for repdate in dates if repdate>=specdate])
    oo=np.array([dd[repdate][specdate]['newCasesLFDOnlyBySpecimenDate'] for repdate in dates if repdate>=specdate])
    cc=np.array([dd[repdate][specdate]['newCasesLFDConfirmedPCRBySpecimenDate'] for repdate in dates if repdate>=specdate])
    print(specdate,"   %5.3f %6d %7.1f"%(1/r[nhist],cc[-1],cc[-1]*(r[nhist]-1)),end='')
    if day-minday>=stabledays:
      c=dccpub[day-minday]# Confirmed by pub date
      d=removedadj2[day-minday]# Denied by pub date
      print("       %5.3f %6d %6d    "%(c/(c+d),c,d),end='')
    else:
      print("           -      -      -    ",end='')
    ccn=cc/nn;oon=oo/nn
    for i in range(min(len(nn)-1,5)):
      ppv=(ccn[-1]-ccn[i])/(oon[i]-oon[-1])
      print(" %5.3f"%ppv,end='')
    print(" ",r[:nhist]*1000,r[1:nhist]/r[0]*100)

"""
py LFDPCRinference.py > LFDPCRoutput
py 7dma.py 2 3 4 5 6 7 < LFDPCRoutput > LFDPCRoutput.7dma
gnuplot
set xdata time;set timefmt "%Y-%m-%d";set format x "%Y-%m-%d"
set xtics rotate by 45 right offset 0.5,0
set terminal pngcairo font "sans,13" size 2560,1280
set bmargin 6;set lmargin 14;set rmargin 15;set tmargin 5
set output "LFDPCR.png"
set xtics "2020-01-06", 604800
set grid xtics lc rgb "#dddddd" lt 1
set title "Estimated probability that a positive LFD test in England that is PCR-retested will be PCR positive"
plot "LFDPCRoutput" u 1:(($3)/(($3)+($4))) w linespoints pt 5 title "PPV by specimen date", "LFDPCRoutput.7dma" u 1:(($3)/(($3)+($4))) w linespoints pt 5 title "PPV by specimen date (7 day moving average)", "LFDPCRoutput" u 1:(($6)/(($6)+($7))) w linespoints pt 5 title "PPV by publication date", "LFDPCRoutput.7dma" u 1:(($6)/(($6)+($7))) w linespoints pt 5 title "PPV by publication date (7 day moving average)"
"""
