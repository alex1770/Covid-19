#!/bin/bash

wget --no-check-certificate https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz -O temp1.csv.gz
if [[ ! -f temp.csv.gz || `diff -q temp1.csv.gz temp.csv.gz` ]]; then
    mv temp1.csv.gz temp.csv.gz
    gunzip -c temp.csv.gz| python sortcog.py > cog_metadata_sorted.csv
    #git add cog_metadata.csv 
    #git commit -m "Update" --date=`date -Iseconds --reference cog_metadata.csv`
fi
