#!/bin/bash

echo
echo HOURLY
date

set -e
python3 getzoemap.py
set +e

python3 getzoenewcases.py
if [ "$?" -eq 1 ]; then
    python3 proczoenewcases.py `date -I`
    cd zoedatapage
    d=`date -I`;x=zoenewcases.London
    cp London.2020-11-05 $x; python3 readincidence.py |grep London >> $x
    cd ..
    python3 proczoenewcases.py $d $x "in London" "per million "
    rsync -a zoenewcases.deconvolve.png zoenewcases.deconvolve.small.png zoenewcases.London.deconvolve.png zoenewcases.London.deconvolve.small.png sonorous@sonorouschocolate.com:public_html/zoe
fi
