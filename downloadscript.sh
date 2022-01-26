#!/bin/bash

set -e

if [ ! -e lastdatadownload ] || [ "`find lastdatadownload -cmin +30`" == "lastdatadownload" ]; then
  for x in confirmed deaths recovered; do
    fn=time_series_covid19_${x}_global.csv
    rm -f $fn
    wget https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/$fn
  done
  #wget https://covid.ourworldindata.org/data/ecdc/full_data.csv -O ecdc.csv
  python3 getworldometerdata.py > worldometer.csv
  #python3 checkdataconsistency.py
  touch lastdatadownload
fi
