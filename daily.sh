#!/bin/bash

echo
echo DAILY
date

set -e

bash regionalcasesbyage.sh

(cd COG-UK; bash update.sh)
(
    cd VOCgrowth
    #python3 uk_var_comp.py BA.1+,BA.1.1+,BA.2+
    #python3 uk_var_comp.py BA.1+,BA.2+
    #python3 uk_var_comp.py BA.1.1+,BA.2+
    #python3 uk_var_comp.py BA.1+,BA.1.1+
    #python3 uk_var_comp.py BA.2,BA.2.1,BA.2.3
    python3 uk_var_comp.py BA.2,BA.2.12.1
    python3 uk_var_comp.py BA.2,BA.4,BA.5
    python3 uk_var_comp.py BA.2,XE
    python3 uk_var_comp.py BA.2,BA.2.18
    python3 uk_var_comp.py BA.2,BA.2.23
)

make
python3 maketrend.py

#. rsynctraffic.sh

#(cd Traffic; python3 parsetraffic.py; python3 maketrafficgraph.py; convert trafficgraph.png -resize '1200x320!' trafficgraph.small.png)

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png recent_cases_growth.png recent_deaths_growth.png'
bigpics=$bigpics' VOCgrowth/UK_BA.2_BA.2.12.1.png VOCgrowth/UK_BA.2_XE.png VOCgrowth/UK_BA.2_BA.2.18.png VOCgrowth/UK_BA.2_BA.2.23.png VOCgrowth/UK_BA.2_BA.4_BA.5.png'
pics=$bigpics
for x in $bigpics; do
    small=${x/.png/.small.png}
    pics="$pics $small"
    convert $x -resize 47% $small
done

rsync -pt worldometer.csv $pics sonorous@sonorouschocolate.com:public_html/covid19/extdata

now=`date -Iminutes`
google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --dump-dom 'https://covid.joinzoe.com/data' > zoedatapage/data.$now 2>> zoedatapage/errors
wget -nv https://covid-assets.joinzoe.com/latest/covid_symptom_study_report.pdf -O zoedatapage/covid_symptom_study_report.$now.pdf 2>> zoedatapage/errors

wget -nv https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data -O VOC/VOCcases.$now 2>> VOC/errors

today=`date -I`
wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=msoa&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDateRollingRate&metric=newCasesBySpecimenDateChange&metric=newCasesBySpecimenDateChangePercentage&format=csv' -O - | gzip > MSOA/msoa_$today.csv.gz
