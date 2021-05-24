#!/bin/bash

echo
echo DAILY
date

set -e

make
python3 maketrend.py

#. rsynctraffic.sh

#(cd Traffic; python3 parsetraffic.py; python3 maketrafficgraph.py; convert trafficgraph.png -resize '1200x320!' trafficgraph.small.png)

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png'
pics=$bigpics
for x in $bigpics; do
    small=${x/.png/.small.png}
    pics="$pics $small"
    convert $x -resize 50% $small
done

rsync -pt worldometer.csv $pics sonorous@sonorouschocolate.com:public_html/covid19/extdata

now=`date -Iminutes`
google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --dump-dom 'https://covid.joinzoe.com/data' > zoedatapage/data.$now 2>> zoedatapage/errors
wget -nv https://covid-assets.joinzoe.com/latest/covid_symptom_study_report.pdf -O zoedatapage/covid_symptom_study_report.$now.pdf 2>> zoedatapage/errors

wget -nv https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data -O VOC/VOCcases.$now 2>> VOC/errors
