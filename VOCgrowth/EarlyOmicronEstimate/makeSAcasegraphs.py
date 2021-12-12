from stuff import *

cases=loadcsv('SAcasecounts.csv')
cases['South Africa']=cases['Total']
del cases['Total']

provinces=['Eastern Cape','Free State','Gauteng','KwaZulu-Natal','Limpopo','Mpumalanga','North West','Northern Cape','Western Cape']

mindate='2021-10-01'
maxdate='3000-12-31'
if len(sys.argv)>1: mindate=sys.argv[1]
if len(sys.argv)>2: maxdate=sys.argv[2]

minday=datetoday(mindate)
maxday=datetoday(maxdate)
day0=datetoday(cases['Date'][0])
minind=max(minday-day0,0)
maxind=min(maxday+1-day0,len(cases['Date']))
for x in cases: cases[x]=cases[x][minind:maxind]
N=len(cases['Date'])
day0=datetoday(cases['Date'][0])
outputdir='output_cases'
os.makedirs(outputdir,exist_ok=True)
dates=[daytodate(day0+d) for d in range(N)]

text=[]
text.append("Data source: https://www.nicd.ac.za/ at "+cases['Date'][-1])

for locname in provinces+['South Africa']:
  cases1=[max(x,0) for x in cases[locname]]
  adjusted=weekdayadj(cases1)
  #adjusted=weekdayadj_slow(cases1)
  
  # Make graphs
  for adj in [0,1]:
    if adj: yy=adjusted
    else: yy=cases1
    
    title='Covid-19 case count in '+locname+' (with weekday adjustment)'*adj
    data=[]
    data.append({
      'title': 'Reported case count'+adj*' (weekday adjusted)',
      'values': list(zip(dates,yy)),
      'with': ('points',1),
      'extra': 'pt 5'
    })
    label='set label at graph 0.4,0.98 "'+'\\n'.join(text)+'"'
    outfn=os.path.join(outputdir,locname.replace(' ','_')+'_cases'+'_adj'*adj+'.png')
    makegraph(title=title, data=data, ylabel='New cases per day', outfn=outfn, extra=[label,'set key left','set logscale y 2'])
