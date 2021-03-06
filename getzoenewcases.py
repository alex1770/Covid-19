import os,sys,requests,json,time

def convdate(x):
  t=time.strptime(x,'%d %B %Y')
  return time.strftime('%Y-%m-%d',t)

tdir='zoedatapage'

def updatenewcases():
  fn=os.path.join(tdir,'zoenewcases')
  last=None
  with open(fn,'r') as fp:
    for x in fp:
      last=x.strip().split()[0]
  r=requests.get('https://covid-assets.joinzoe.com/latest/incidence.json')
  if r.status_code!=200: raise ConnectionError("Error code %d from new cases api"%r.status_code)
  js=r.json()
  dt=convdate(js['uk_incidence_updated_on'])
  if last==None or dt>last:
    n=int(js['uk_incidence'].replace(',',''))
    with open(fn,'a') as fp:
      print(dt,"%7d"%n,file=fp)
    print("Wrote new cases:",n)
    return True
  return False
  
if __name__=="__main__":
  sys.exit(not updatenewcases())
