#!/bin/bash

if [ $# -ne 2 ] ; then
    echo "Usage: dir_name file_name_prefix"
    echo "Output: paths of files that has the same job id."
    exit
fi

DIR_NAME=$1
PREFIX=$2

TOTAL_NUMBER_OF_FILES=$( ls $DIR_NAME | wc -l )

echo "Total number of files in the directory: $TOTAL_NUMBER_OF_FILES"

FILE_LIST=$( ls $DIR_NAME | sed "s/\(${PREFIX}_\)\([0-9]*\)\(.*\)/\2/" | sort | uniq -d )

if [ "x$FILE_LIST" = "x" ] ; then
	echo "There are no jobs that have more than one output file."
else
	echo "IDs of the jobs that have more than one output file:" $FILE_LIST
	echo "$FILE_LIST" | sed "s/\(.*\)/${PREFIX}_\1_\*/" | xargs -n 1 find $DIR_NAME -name | xargs ls -al
fi
