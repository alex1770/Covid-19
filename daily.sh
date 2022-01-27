#!/bin/bash

echo
echo DAILY
date

set -e

bash regionalcasesbyage.sh

(cd COG-UK; bash update.sh)
(cd VOCgrowth; python3 uk_var_comp.py; python3 uk_var_comp.py BA.1,BA.2; python3 uk_var_comp.py BA.1.1,BA.2)

make
python3 maketrend.py

#. rsynctraffic.sh

#(cd Traffic; python3 parsetraffic.py; python3 maketrafficgraph.py; convert trafficgraph.png -resize '1200x320!' trafficgraph.small.png)

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png recent_cases_growth.png recent_deaths_growth.png VOCgrowth/UK_BA.1_BA.1.1_BA.2.png VOCgrowth/UK_BA.1_BA.2.png VOCgrowth/UK_BA.1.1_BA.2.png'
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

today=`date -I`
wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=msoa&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDateRollingRate&metric=newCasesBySpecimenDateChange&metric=newCasesBySpecimenDateChangePercentage&format=csv' -O - | gzip > MSOA/msoa_$today.csv.gz
