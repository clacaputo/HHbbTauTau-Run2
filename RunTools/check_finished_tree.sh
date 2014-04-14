#!/bin/bash

if [ $# -ne 3 ] ; then
    echo "Usage: file_list_path job_result new_file_list_path"
    exit
fi

FILE_LIST_PATH=$1
FILE_JOB_RESULT=$2
NEW_FILE_LIST_PATH=$3

if [ ! -d "$FILE_LIST_PATH" ] ; then
        echo "ERROR: file list path '$FILE_LIST_PATH' does not exist."
        exit
fi

JOBS=$( find $FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" | sed "s/\.txt//" | sort )

if [ "x$JOBS" = "x" ] ; then
        echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
        exit
fi

SUCCESSFULL_JOBS=$( cat $FILE_JOB_RESULT | sed -n "s/\(^0 \)\([^ ]*\)\(.*\)/\2/p" | sort )

echo -e "$JOBS"\\n"$SUCCESSFULL_JOBS" | sort | uniq -u | xargs -n 1 printf "$FILE_LIST_PATH/%b.txt $NEW_FILE_LIST_PATH/\n" \
     | xargs -n 2 cp 
