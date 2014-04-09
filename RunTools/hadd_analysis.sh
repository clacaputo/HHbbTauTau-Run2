#!/bin/bash

if [ $# -ne 2 ] ; then
    echo "Usage: analysis_path output_dir"
    exit
fi

ANALYSIS_PATH=$1
cd $ANALYSIS_PATH

FOLDER_PATH="anaMuTau anaETau anaTauTau"
OUTPUT_DIR=$2

if [ -d "$OUTPUT_DIR" ] ; then
	echo "ERROR: output dir '$OUTPUT_DIR' already exists."
	exit
fi

mkdir -p $OUTPUT_DIR

for FOLDER in $FOLDER_PATH ; do
    SUB_FOLDERS=$( find $FOLDER -maxdepth 1 -type d -printf "%f\n" )
    for SUB_FOLDER in $SUB_FOLDERS ; do
		if [ $SUB_FOLDER = $FOLDER ] ; then
			continue
		fi
		hadd $OUTPUT_DIR/${FOLDER}_${SUB_FOLDER}.root $FOLDER/$SUB_FOLDER/*.root
	done
done

tar cjvf ${OUTPUT_DIR}.tar.bz2 $OUTPUT_DIR
