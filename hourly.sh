#!/bin/bash

echo
echo HOURLY
date

set -e
python3 getzoemap.py
set +e

python3 getzoenewcases.py
if [ "$?" -eq 1 ]; then
    python3 proczoenewcases.py
    rsync -a zoenewcasesdeconvolve.png zoenewcasesdeconvolve.small.png sonorous@sonorouschocolate.com:public_html/zoe
fi
