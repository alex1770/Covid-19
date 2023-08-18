import requests,bs4
import sys

# Make computer-readable alias list from cov-lineages file https://cov-lineages.org/lineage_list.html
# Superseded by https://github.com/cov-lineages/pango-designation/raw/master/pango_designation/alias_key.json (includes recombinants)
#               (from https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json)
#               or https://cov-spectrum.org/api/v2/resource/pango-lineage-alias
# 
# (Note that https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt has further details)
page=requests.get("https://cov-lineages.org/lineage_list.html")
soup=bs4.BeautifulSoup(page.content,features="lxml")

d={}

for x in soup.find_all('tr'):
  lin=None;alias=None
  for y in x.find_all('td'):
    if lin==None: lin=y.text.strip()
    if y.text[:5]=="Alias":
      alias=y.text.split()[2].rstrip(',')
      if alias=="Alias": alias=y.text.split()[4].rstrip(',')# Workaround typo in page source
  if lin!=None and alias!=None:
    p=lin.find('.')
    assert p>=0
    if lin[p:]!=alias[-(len(lin)-p):]:
      print("Inconsistent alias",lin,alias,file=sys.stderr)
      continue
    alias=alias[:-(len(lin)-p)];lin=lin[:p]
    if lin in d: assert d[lin]==alias
    else: d[lin]=alias

# Sort in reverse length order so that an earlier alias can't be a prefix of a later one
l=list(d)
l.sort(key=lambda x:-len(d[x]))

with open("variantaliases.py","w") as fp:
  print("# Automatically generated using makealiaslist.py",file=fp)
  print("aliases=[",file=fp)
  for x in l:
    print('("%s","%s"),'%(x,d[x]),file=fp)
  print("]",file=fp)
  
