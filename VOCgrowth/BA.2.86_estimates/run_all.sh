#!/bin/bash

# Get directory path
DIR="./counts"

# Check if directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory $DIR does not exist."
    exit 1
fi

# Collect all the files in the directory
files=("$DIR"/*)

# Execute the python script with all the files as arguments
python binomialfit.py "${files[@]}"
