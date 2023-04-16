#!/bin/bash

echo
echo DAILY
date

set -e

#bash regionalcasesbyage.sh

(cd COG-UK; bash update.sh)
(
    cd VOCgrowth
    python3 uk_var_comp.py -f 2023-02-01 -p -b -l 'XBB.1.5*,XBB.1.9*,XBB.1.16*' -F30
    python3 uk_var_comp.py -f 2023-02-01 -l 'BQ.1,BQ.1.1,BQ.1.8,BQ.1.1.2,BQ.1.1.8,CH.1.1,CH.1.1.1,XBF,DN.1.1,EG.1,EM.1,DU.1,XBB.1,*,DV.1*,XBB.*,XBB.1.*,XBB.1.5*,XBB.1.9.1*,XBB.1.9.2*,XBB.1.16*,BA.5*,BA.5.2*,BF.*,BN.1.*,BQ.1.*,BQ.1.1*' -F30
)

make
python3 maketrend.py

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png recent_cases_growth.png recent_deaths_growth.png'
bigpics=$bigpics' VOCgrowth/UK_XBB.1.5*_XBB.1.9*_XBB.1.16*.png'
bigpics=$bigpics' VOCgrowth/UK_BQ.1_BQ.1.1_BQ.1.8_BQ.1.1.2_BQ.1.1.8_CH.1.1_CH.1.1.1_XBF_DN.1.1_EG.1_EM.1_DU.1_XBB.1_*_DV.1*_XBB.*_XBB.1.*_XBB.1.5*_XBB.1.9.1*_XBB.1.9.2*_XBB.1.16*_BA.5*_BA.5.2*_BF.*_BN.1.*_BQ.1.*_BQ.1.1*.variantpressure.png'
bigpics=$bigpics' VOCgrowth/UK_BQ.1_BQ.1.1_BQ.1.8_BQ.1.1.2_BQ.1.1.8_CH.1.1_CH.1.1.1_XBF_DN.1.1_EG.1_EM.1_DU.1_XBB.1_*_DV.1*_XBB.*_XBB.1.*_XBB.1.5*_XBB.1.9.1*_XBB.1.9.2*_XBB.1.16*_BA.5*_BA.5.2*_BF.*_BN.1.*_BQ.1.*_BQ.1.1*.growthproj.png'

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
