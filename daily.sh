#!/bin/bash

echo
echo DAILY
date

set -e

#bash regionalcasesbyage.sh

(cd COG-UK; bash update.sh)
(
    cd VOCgrowth
    python3 uk_var_comp.py -f 2023-05-01 -p -b -l 'XBB.1.16,XBB.1.16.11,EG.5.1*' -F30
    python3 -i uk_var_comp.py -f 2023-05-01 -l 'CH.1.1.1,EG.1,EG.5.1,EG.5.1.1,FU.1,GE.1,XBB.1.16,XBB.1.16.11,XBB.2.3,XBB.2.3.2,XBB.2.3.3,*,XBB.*,XBB.1.*,XBB.1.5*,XBB.1.9.1*,XBB.1.9.2*,XBB.1.16*,BQ.1.1*,EG.*,FL.*,FE.*,FY.*' -F30
)

make
python3 maketrend.py

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png recent_cases_growth.png recent_deaths_growth.png'
bigpics=$bigpics' VOCgrowth/UK_XBB.1.16_XBB.1.16.11_EG.5.1*.png'
bigpics=$bigpics' VOCgrowth/UK_CH.1.1.1_EG.1_EG.5.1_EG.5.1.1_FU.1_GE.1_XBB.1.16_XBB.1.16.11_XBB.2.3_XBB.2.3.2_XBB.2.3.3_*_XBB.*_XBB.1.*_XBB.1.5*_XBB.1.9.1*_XBB.1.9.2*_XBB.1.16*_BQ.1.1*_EG.*_FL.*_FE.*_FY.*.variantpressure.png'
bigpics=$bigpics' VOCgrowth/UK_CH.1.1.1_EG.1_EG.5.1_EG.5.1.1_FU.1_GE.1_XBB.1.16_XBB.1.16.11_XBB.2.3_XBB.2.3.2_XBB.2.3.3_*_XBB.*_XBB.1.*_XBB.1.5*_XBB.1.9.1*_XBB.1.9.2*_XBB.1.16*_BQ.1.1*_EG.*_FL.*_FE.*_FY.*.growthproj.png'

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
