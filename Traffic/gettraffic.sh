#!/bin/bash

fn=London.`date -Iminutes`.png
/usr/bin/google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --screenshot=$fn --window-size=1600,1200 'https://www.google.co.uk/maps/@51.5121089,-0.1476081,11z/data=!5m1!1e1?hl=en' 2>> errors
