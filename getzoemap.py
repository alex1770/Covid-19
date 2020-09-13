import os
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

opts = Options()
opts.add_argument("--headless")

driver = webdriver.Chrome(chrome_options=opts)
#driver = webdriver.Chrome('/usr/bin/chromedriver',chrome_options=opts)

driver.get('file://'+os.path.join(os.getcwd(),'getzoemap.html'))

WebDriverWait(driver, 100).until(EC.presence_of_element_located((By.ID, 'nowfinished')))

res=driver.find_element_by_id('xyzzy')
a=res.text

driver.quit()

print(a)
