import os,json,sys
from requests import get

def get_data(req):
  url='https://api.coronavirus.data.gov.uk/v2/data?'
  response = get(url+req, timeout=10)
  if not response.ok:
    raise RuntimeError(f'Request failed: { response.text }')
  data=response.json()['body'][::-1]
  for d in data:
    for x in d:
      if d[x]==None: d[x]=0
  return data

dirname='apidata_lfd'

req='areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDate&metric=newCasesLFDConfirmedPCRBySpecimenDate&metric=newCasesLFDOnlyBySpecimenDate&metric=newLFDTests&format=json'

data=get_data(req)

updatedate=data[-1]['date']

fn=os.path.join(dirname,updatedate)
if len(sys.argv)==1 and os.path.isfile(fn): sys.exit(1)# Exit signalling no update needs to be done
os.makedirs(dirname, exist_ok=True)
with open(fn,'w') as fp:
  json.dump(data,fp,indent=2)
