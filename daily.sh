#!/bin/bash

echo
echo DAILY
date

set -e

#bash regionalcasesbyage.sh

(cd COG-UK; bash update.sh)
(
    cd VOCgrowth
    #python3 uk_var_comp.py BA.2,BA.2.12.1
    #python3 uk_var_comp.py BA.2,BA.4,BA.5,BA.5.1
    #python3 uk_var_comp.py -f 2022-05-01 -p -l 'BA.2*,BA.2.12.1,BA.4*,BA.5*'
    #python3 uk_var_comp.py -f 2022-06-01 -b -l 'BA.5.1,BA.4,BE.1,BA.5.2,BA.5.2.1'
    #python3 uk_var_comp.py -f 2022-05-01 -b -p -l 'BA.5.1,BA.2.12.1,BA.4'
    python3 uk_var_comp.py -f 2022-06-01 -b -p -l 'BA.5.1,BA.4.6,BA.5.2,BA.5.2.1,BF.5'
    python3 uk_var_comp.py -f 2022-06-01 -b -p -l 'BA.5*,BA.2.75'
    python3 uk_var_comp.py -f 2022-06-01 -b -p -l 'BA.5*,BA.2.75,BA.2.76'
    python3 uk_var_comp.py -f 2022-06-01 -b -l 'BA.5.1,BA.4.6,BA.5.2,BA.5.2.1,BF.5,BA.2*,BA.4*,BA.5*' -f1 50 -f2 50
)

make
python3 maketrend.py

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png recent_cases_growth.png recent_deaths_growth.png'
bigpics=$bigpics' VOCgrowth/UK_BA.5.1_BA.4.6_BA.5.2_BA.5.2.1_BF.5.png VOCgrowth/UK_BA.5*_BA.2.75.png VOCgrowth/UK_BA.5*_BA.2.75_BA.2.76.png VOCgrowth/UK_BA.5.1_BA.4.6_BA.5.2_BA.5.2.1_BF.5_BA.2*_BA.4*_BA.5*.png VOCgrowth/UK_BA.5.1_BA.4.6_BA.5.2_BA.5.2.1_BF.5_BA.2*_BA.4*_BA.5*.growthproj.png'
set -o noglob
pics=$bigpics
for x in $bigpics; do
    small=${x/.png/.small.png}
    pics="$pics $small"
    cat $x | convert - -resize 47% - > $small
done
pics="$pics VOCgrowth/UK_BA.5*_BA.2.75_BA.2.76"
set +o noglob

rsync -pt worldometer.csv $pics sonorous@sonorouschocolate.com:public_html/covid19/extdata

now=`date -Iminutes`
google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --dump-dom 'https://covid.joinzoe.com/data' > zoedatapage/data.$now 2>> zoedatapage/errors
wget -nv https://covid-assets.joinzoe.com/latest/covid_symptom_study_report.pdf -O zoedatapage/covid_symptom_study_report.$now.pdf 2>> zoedatapage/errors

# wget -nv https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data -O VOC/VOCcases.$now 2>> VOC/errors

today=`date -I`
wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=msoa&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDateRollingRate&metric=newCasesBySpecimenDateChange&metric=newCasesBySpecimenDateChangePercentage&format=csv' -O - | gzip > MSOA/msoa_$today.csv.gz
