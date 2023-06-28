#!/bin/bash

fn=London.`date -Iminutes`.png
python3 seleniumscreenshot.py "$fn"
