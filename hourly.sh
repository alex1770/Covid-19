#!/bin/bash

echo
echo HOURLY
date

set -e
python3 getzoemap.py
set +e

# Trigger (lightweight) for new incidence (swab-based) data
python3 getzoenewcases.py
if [ "$?" -eq 0 ]; then
    now=`date -I`
    python3 proczoenewcases.py $now
    cd zoedatapage
    google-chrome-stable --headless --disable-gpu --disable-features=NetworkService --dump-dom 'https://covid-assets.joinzoe.com/latest/incidence_map.html?v=1.1' > incidence.$now 2>> errors
    x=zoenewcases.London
    cp London.2020-11-05 $x; python3 readincidence.py |grep London >> $x
    cd ..
    python3 proczoenewcases.py $now $x "in London" "per million "
    rsync -a zoenewcases.deconvolve.png zoenewcases.deconvolve.small.png zoenewcases.London.deconvolve.png zoenewcases.London.deconvolve.small.png sonorous@sonorouschocolate.com:public_html/zoe
fi
