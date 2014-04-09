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
        mkdir -p $OUTPUT_DIR/$FOLDER
        if [ $SUB_FOLDER = "Radion" ] ; then
            cp $FOLDER/$SUB_FOLDER/*.root $OUTPUT_DIR/$FOLDER/
        else
            hadd $OUTPUT_DIR/$FOLDER/${SUB_FOLDER}.root $FOLDER/$SUB_FOLDER/*.root
        fi
	done
done

tar cjvf ${OUTPUT_DIR}.tar.bz2 $OUTPUT_DIR
