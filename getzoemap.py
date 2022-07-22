#/usr/bin/env python3

import sys,os,ast,json
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

#from proczoemap import processdata
#from proczoemap_reg import processdata_reg
from subprocess import Popen

floatkeys=["long", "lat", "st_areasha", "st_lengths", "corrected_covid_positive", "cases", "cases_pm", "percentage", "discrete_percentage", "respondent", "predicted_covid_positive_count"]
intkeys=["cartodb_id", "objectid", "bng_e", "bng_n", "population", "discrete_cases_pm"]

opts=Options()
opts.add_argument("--headless")
driver=webdriver.Chrome(options=opts)
#driver=webdriver.Chrome('/usr/bin/chromedriver',options=opts)
driver.get('file://'+os.path.join(os.getcwd(),'getzoemap.html'))
WebDriverWait(driver,100).until(EC.presence_of_element_located((By.ID, 'nowfinished')))
res=driver.find_element_by_id('xyzzy').text.split('\n')
driver.quit()

# Convert format HH:MM:SS DD-MM-YYYY to YYYY-MM-DD
def convdate(ds):
  return ds[-4:]+'-'+ds[-7:-5]+'-'+ds[-10:-8]

tdir='zoemapdata'
date=convdate(ast.literal_eval(res[0])['data_status'])
fn=os.path.join(tdir,date)
if os.path.isfile(fn): sys.exit(1)
d={}
for r in res:
  e=ast.literal_eval(r)
  for k in list(e):
    if k in floatkeys:
      try:
        e[k]=float(e[k])
      except:
        del e[k]
    elif k in intkeys:
      try:
        e[k]=int(e[k])
      except:
        del e[k]
  d[e["utla19nm"]]=e
with open(fn,'w') as fp:
  json.dump(d,fp,indent=2)
#processdata(tdir)
#processdata_reg(tdir)
#Popen("rsync -a zoeselected.csv zoeselected.png zoeregions.csv zoeregions.png sonorous@sonorouschocolate.com:public_html/zoe",shell=True).wait()
sys.exit(0)
