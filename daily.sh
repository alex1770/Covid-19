#!/bin/bash

echo
echo DAILY
date

set -e

#bash regionalcasesbyage.sh

(cd COG-UK; bash update.sh)
(
    cd VOCgrowth
    python3 uk_var_comp.py -f 2023-07-01 -p -b -l 'EG.5.1*,XBB.1.16*,BA.2.86*' -F30
    python3 uk_var_comp.py -f 2023-07-01 -l 'EG.1,EG.5.1,EG.5.1.1,EG.5.1.3,EG.5.1.4,EG.5.1.6,GE.1,XBB.1.16,XBB.1.16.11,XBB.1.16.15,XBB.1.16.6,XBB.2.3,XBB.2.3.11,DV.7.1,FL.1.5.1,GK.2,*,XBB.*,XBB.1.*,XBB.1.5*,XBB.1.9.*,XBB.1.16*,EG.*,FL.*,FY.*,BA.2.86*' -F30
)

make
python3 maketrend.py

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png recent_cases_growth.png recent_deaths_growth.png'
bigpics=$bigpics' VOCgrowth/UK_EG.5.1*_XBB.1.16*_BA.2.86*.png'
bigpics=$bigpics' VOCgrowth/UK_EG.1_EG.5.1_EG.5.1.1_EG.5.1.3_EG.5.1.4_EG.5.1.6_GE.1_XBB.1.16_XBB.1.16.11_XBB.1.16.15_XBB.1.16.6_XBB.2.3_XBB.2.3.11_DV.7.1_FL.1.5.1_GK.2_*_XBB.*_XBB.1.*_XBB.1.5*_XBB.1.9.*_XBB.1.16*_EG.*_FL.*_FY.*_BA.2.86*.variantpressure.png'
bigpics=$bigpics' VOCgrowth/UK_EG.1_EG.5.1_EG.5.1.1_EG.5.1.3_EG.5.1.4_EG.5.1.6_GE.1_XBB.1.16_XBB.1.16.11_XBB.1.16.15_XBB.1.16.6_XBB.2.3_XBB.2.3.11_DV.7.1_FL.1.5.1_GK.2_*_XBB.*_XBB.1.*_XBB.1.5*_XBB.1.9.*_XBB.1.16*_EG.*_FL.*_FY.*_BA.2.86*.growthproj.png'

set -o noglob
pics=$bigpics
for x in $bigpics; do
    small=${x/.png/.small.png}
    pics="$pics $small"
    cat $x | convert - -resize 47% - > $small
done
set +o noglob

rsync -pt worldometer.csv $pics sonorous@sonorouschocolate.com:public_html/covid19/extdata

now=`date -Iminutes`
google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --dump-dom 'https://covid.joinzoe.com/data' > zoedatapage/data.$now 2>> zoedatapage/errors
wget -nv https://covid-assets.joinzoe.com/latest/covid_symptom_study_report.pdf -O zoedatapage/covid_symptom_study_report.$now.pdf 2>> zoedatapage/errors

# wget -nv https://www.gov.uk/government/publications/covid-19-variants-genomically-confirmed-case-numbers/variants-distribution-of-cases-data -O VOC/VOCcases.$now 2>> VOC/errors

today=`date -I`
wget 'https://api.coronavirus.data.gov.uk/v2/data?areaType=msoa&metric=newCasesBySpecimenDateRollingSum&metric=newCasesBySpecimenDateRollingRate&metric=newCasesBySpecimenDateChange&metric=newCasesBySpecimenDateChangePercentage&format=csv' -O - | gzip > MSOA/msoa_$today.csv.gz
