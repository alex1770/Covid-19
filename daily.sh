#!/bin/bash

echo
echo DAILY
date

set -e

make
python3 maketrend.py

python3 hosp+casesbyage.py

#. rsynctraffic.sh

#(cd Traffic; python3 parsetraffic.py; python3 maketrafficgraph.py; convert trafficgraph.png -resize '1200x320!' trafficgraph.small.png)

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png hospitaladmissionsbyage-abs.png admissionandcaseageratios.png'
pics=$bigpics
for x in $bigpics; do pics="$pics ${x/.png/.small.png}"; done

for x in $bigpics; do
  convert $x -resize 50% ${x/.png/.small.png}
done

rsync -pt worldometer.csv $pics sonorous@sonorouschocolate.com:public_html/covid19/extdata

now=`date -Iminutes`
google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --dump-dom 'https://covid.joinzoe.com/data' > zoedatapage/data.$now 2>> zoedatapage/errors
wget -nv https://covid-assets.joinzoe.com/latest/covid_symptom_study_report.pdf -O zoedatapage/covid_symptom_study_report.$now.pdf 2>> zoedatapage/errors
