#!/bin/bash

dir=casesbyageregion
mkdir -p $dir
skip=1
(
    for loc in England North_East North_West Yorkshire_and_The_Humber East_Midlands West_Midlands East_of_England London South_East South_West; do
        (
            python3 casesbyage.py -l $loc -s $skip -o $dir/$loc.png
            convert $dir/$loc.png -resize 20% $dir/$loc.small.png
        ) &
    done
    wait
    rsync -a $dir sonorous@sonorouschocolate.com:public_html/covid19/extdata/
)
