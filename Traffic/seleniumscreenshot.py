import sys,time

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Set Chrome to run headlessly
chrome_options = Options()
chrome_options.add_argument("--headless")

# Initialize the driver
driver = webdriver.Chrome(options=chrome_options)

# Set screen size
driver.set_window_size(1600, 1200)  # Set the width and height

# Navigate to the site
driver.get('https://www.google.co.uk/maps/@51.5121089,-0.1476081,11z/data=!5m1!1e1?hl=en')

# wait for the consent form to load
wait = WebDriverWait(driver, 10)
accept_button = wait.until(EC.element_to_be_clickable((By.XPATH, "//*[text()='Accept all']")))

# click on the "accept" button
accept_button.click()

time.sleep(10)# Probably 1 second is plenty, but who cares?
#WebDriverWait(driver, 10).until(lambda d: d.execute_script('return jQuery.active == 0'))

# Take a screenshot
driver.save_screenshot(sys.argv[1])

# Quit the driver
driver.quit()
