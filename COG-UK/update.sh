#!/bin/bash

wget https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv.gz -O temp.csv.gz
if [[ ! -f cog_metadata.csv.gz || `diff -q temp.csv.gz cog_metadata.csv.gz` ]]; then
    mv temp.csv.gz cog_metadata.csv.gz
    gunzip -c cog_metadata.csv.gz| python sortcog.py > cog_metadata_sorted.csv
    #git add cog_metadata.csv 
    #git commit -m "Update" --date=`date -Iseconds --reference cog_metadata.csv`
fi
