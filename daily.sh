#!/bin/bash

echo
echo DAILY
date

set -e

make
python3 maketrend.py

#. rsynctraffic.sh

#(cd Traffic; python3 parsetraffic.py; python3 maketrafficgraph.py; convert trafficgraph.png -resize '1200x320!' trafficgraph.small.png)

for x in trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png; do
  convert $x -resize 50% ${x/.png/.small.png}
done

rsync -pt worldometer.csv trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendthr_cases.small.png trendthr_deaths.small.png trendsimple_cases.small.png trendsimple_deaths.small.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png trendsimple_cases_zoom.small.png trendsimple_deaths_zoom.small.png recent_cases.small.png recent_cases.png recent_deaths.png recent_deaths.small.png sonorous@sonorouschocolate.com:public_html/covid19/extdata

now=`date -Iminutes`
google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --dump-dom 'https://covid.joinzoe.com/data' > zoedatapage/data.$now 2>> zoedatapage/errors
wget -nv https://covid-assets.joinzoe.com/latest/covid_symptom_study_report.pdf -O zoedatapage/covid_symptom_study_report.$now.pdf 2>> zoedatapage/errors
