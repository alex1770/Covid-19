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
    #python3 uk_var_comp.py -f 2022-06-01 -b -p -l 'BA.5.1,BA.4.6,BA.5.2,BA.5.2.1,BF.5'
    #python3 uk_var_comp.py -f 2022-06-01 -b -p -l 'BA.5*,BA.2.75,BA.2.75.1'
    #python3 uk_var_comp.py -f 2022-05-01 -b -l 'BA.5.1,BA.4.1,BA.4.6,BA.5.1+S:R346T,BA.5.1.5,BA.5.5,BA.5.2,BA.5.2.1,BA.5.2.6,BA.5.2.7,BA.5.3.3,BE.1,BE.1.1,BF.3,BF.5,BF.6,BF.7,BF.10,BF.11,BA.2.12.1,BA.2.75,BA.2.75.1,BA.2*,BA.4*,BA.5*' -F 30
    #python3 uk_var_comp.py -f 2022-05-01 -l 'BA.5.1,BA.2.75.1,BA.2.75.5,BA.4.6,BA.5.1.12,BA.5.1+S:R346T,BA.5.2,BA.5.2.1,BA.5.2.6,BA.5.2.7,BE.1,BE.1.1,BF.11,BF.3,BF.5,BF.7,BA.2.12.1,BA.1*,BA.2*,BA.4*,BA.5*' -F30
    #python3 uk_var_comp.py -f 2022-07-01 -b -l 'BA.5.2,BF.7,BA.2.75.1' -F30
    #python3 uk_var_comp.py -f 2022-07-01 -b -l 'BA.5.2,BA.2.75.2,BA.2.75.5' -F30
    #python3 uk_var_comp.py -f 2022-07-01 -b -l 'BA.5.2,BF.7,BA.2.75.1,BA.2.75.2,BA.2.75.5' -F30
    #python3 uk_var_comp.py -f 2022-05-01 -l 'BA.5.1,BA.2.75.1,BA.2.75.2,BA.2.75.5,BA.4.6,BA.5.1.12,BA.5.1+S:R346T,BA.5.2,BA.5.2.1,BA.5.2.6,BA.5.2.7,BE.1,BE.1.1,BF.11,BF.3,BF.5,BF.7,BA.2.12.1,BA.1*,BA.2*,BA.4*,BA.5*' -F30
    #python3 uk_var_comp.py -f 2022-05-01 -l 'BA.5.1,BA.2.75.2,BA.2.75.5,BN.1,BA.4.6,BA.5.1.12,BA.5.2,BA.5.2.1,BA.5.2.6,BA.5.2.7,BF.7,BM.1.1,BA.2.12.1,BQ.1,BQ.1.1,BA.1*,BA.2*,BA.4*,BA.5*,BA.2.75.*' -F30
    #python3 uk_var_comp.py -f 2022-09-01 -p -b -l 'BA.5.2,BQ.1,BQ.1.1,XBB.1' -F30
    #python3 uk_var_comp.py -f 2022-09-01 -l 'BA.5.1,BA.2.75.2,BA.2.75.5,BN.1,BA.4.6,BA.5.1.12,BA.5.2,BA.5.2.1,BA.5.2.6,BA.5.2.7,BF.7,BF.11,BM.1.1,BQ.1,BQ.1.1,XBB,XBB.1,BA.*,BA.2.75.*,BA.4*,BA.5*' -F30
    python3 uk_var_comp.py -f 2022-12-01 -p -b -l 'BQ.1.1*,CH.1.1.1,XBB.1.5' -F30
    python3 uk_var_comp.py -f 2022-11-15 -l 'BA.2.75.2,BA.4.6,BA.5.2.6,BA.5.2.13,BA.5.2.35,BF.5,BF.7,BF.11,BQ.1,BQ.1.8,BQ.1.1.8,BE.1.1,CH.1.1,CH.1.1.1,XBB,XBB.1,*,XBB.*,XBB.1.*,BA.2.75.*,BA.4*,BA.5*,BF.*,BN.1.*,BQ.1.*,BQ.1.1*' -F30
)

make
python3 maketrend.py

bigpics='trendthr_cases.png trendthr_deaths.png trendsimple_cases.png trendsimple_deaths.png trendsimple_cases_zoom.png trendsimple_deaths_zoom.png recent_cases.png recent_deaths.png recent_cases_growth.png recent_deaths_growth.png'
bigpics=$bigpics' VOCgrowth/UK_BQ.1.1*_CH.1.1.1_XBB.1.5.png'
bigpics=$bigpics' VOCgrowth/UK_BA.2.75.2_BA.4.6_BA.5.2.6_BA.5.2.13_BA.5.2.35_BF.5_BF.7_BF.11_BQ.1_BQ.1.8_BQ.1.1.8_BE.1.1_CH.1.1_CH.1.1.1_XBB_XBB.1_*_XBB.*_XBB.1.*_BA.2.75.*_BA.4*_BA.5*_BF.*_BN.1.*_BQ.1.*_BQ.1.1*.variantpressure.png'
bigpics=$bigpics' VOCgrowth/UK_BA.2.75.2_BA.4.6_BA.5.2.6_BA.5.2.13_BA.5.2.35_BF.5_BF.7_BF.11_BQ.1_BQ.1.8_BQ.1.1.8_BE.1.1_CH.1.1_CH.1.1.1_XBB_XBB.1_*_XBB.*_XBB.1.*_BA.2.75.*_BA.4*_BA.5*_BF.*_BN.1.*_BQ.1.*_BQ.1.1*.growthproj.png'

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
