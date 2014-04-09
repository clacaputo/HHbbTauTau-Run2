#!/bin/bash

if [ $# -ne 1 ] ; then
    echo "Usage: analysis_path"
    exit
fi

ANALYSIS_PATH=$1
cd $ANALYSIS_PATH

FOLDER_PATH="anaMuTau" "anaETau" "anaTauTau"
OUTPUT_DIR="anaResult"

mkdir -p $OUTPUT_DIR

for FOLDER in $FOLDER_PATH ; do
    SUB_FOLDER=$( find $FOLDER -maxdepth 1 -type d -printf "%f\n" )
    echo $SUB_FOLDER
done
