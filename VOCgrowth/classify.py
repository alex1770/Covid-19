import datetime,os,json

aliasfile="alias_key.json"
aliasmtime=datetime.datetime.utcfromtimestamp(os.path.getmtime(aliasfile))#.strftime('%Y-%m-%d-%H-%M-%S')
now=datetime.datetime.utcnow()

if (now-aliasmtime).seconds>=3600:
  import requests
  page=requests.get("https://github.com/cov-lineages/pango-designation/raw/master/pango_designation/"+aliasfile)
  with open(aliasfile,'w') as fp:
    fp.write(page.text)
  rawaliases=page.json()
else:
  with open(aliasfile) as fp:
    rawaliases=json.load(fp)

# 'A' and 'B' shouldn't alias to empty string. Remove recombinants for the time being.
aliases=[(c,e) for (c,e) in rawaliases.items() if e!="" and type(e)==str]

# Sort in decreasing order of size, so that longest substitution possible is made in contractlin.
aliases.sort(key=lambda x:-len(x[1]))

ecache={}
def expandlin(lin):
  if lin in ecache: return ecache[lin]
  for (short,long) in aliases:
    s=len(short)
    if lin[:s+1]==short+".": ecache[lin]=long+lin[s:];return ecache[lin]
  ecache[lin]=lin
  return lin

ccache={}
def contractlin(lin):
  if lin in ccache: return ccache[lin]
  lin=expandlin(lin)
  for (short,long) in aliases:
    l=len(long)
    if lin[:l+1]==long+".": ccache[lin]=short+lin[l:];return ccache[lin]
  ccache[lin]=lin
  return lin

# Manual classifications are largely superseded by better (Usher) COG classifications and cov-spectrum/Lapis queries for GISAID

# Return classification cl1 if (cl0 is Unassigned) or (cl0 is a prefix of cl1) or (cl1 is a recombinant and cl0 isn't), otherwise return cl0
# So cl0 takes priority if there is a disagreement
def join(cl0,cl1):
  if cl1=="Other": cl1="Unassigned"
  if cl0=="Unassigned": return cl1
  if cl1[0]=="X" and cl0[0]!="X": return cl1
  el0=expandlin(cl0)
  el1=expandlin(cl1)
  if el1[:len(el0)]==el0: return cl1
  return cl0

def classify(mutations,cl="",gisaid=False):
  if gisaid:
    cl=join(cl,classifygisaid.treeclassify(mutations))
    cl=join(cl,manualclassifygisaid(cl,mutations))
  else:
    cl=join(cl,classifycog.treeclassify(mutations))
    cl=join(cl,manualclassifycog(cl,mutations))
  return cl
