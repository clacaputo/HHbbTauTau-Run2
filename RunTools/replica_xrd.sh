#!/bin/bash

if [ $# -lt 2 || $# -gt 3 ] ; then
    echo "Usage: file_list_path output_path [max_n_files]"
    exit
fi

FILE_LIST_PATH=$1
OUTPUT_PATH=$2

if [ $# -eq 3 ] ; then
    MAX_N_FILES=$3
else
    MAX_N_FILES=0
fi

SERVER="root://xrootd.ba.infn.it/"

if [ ! -d $FILE_LIST_PATH ] ; then
    echo "ERROR: file list path '$FILE_LIST_PATH' does not exists."
    exit
fi

if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: output path '$OUTPUT_PATH' does not exist."
    exit
fi

JOB_FILES=$( find $FILE_LIST_PATH -maxdepth 1 -name "*.txt" )
N_JOB_FILES=$( echo "$JOB_FILES" | wc -l )
if [ $N_JOB_FILES -eq 0 ] ; then
    echo "ERROR: no job files found in '$FILE_LIST_PATH'"
    exit
fi

n=0

for JOB_FILE in $JOB_FILES ; do
    FILES_TO_TRANSFER=$( cat $JOB_FILE | sed '/^\s*$/d' )
    N_FILES_TO_TRANSFER=$( echo "$FILES_TO_TRANSFER" | sed '/^\s*$/d' | wc -l )
    if [ $N_FILES_TO_TRANSFER -eq 0 ] ; then
        echo "ERROR: input job file '$JOB_FILE' is empty."
        exit
    fi
    for FILE in $FILES_TO_TRANSFER ; do
        echo "Starting transfering file '$FILE'..."
        while [ 1 ] ; do
            xrdcp -force "${SERVER}${FILE}" $OUTPUT_PATH
            RESULT=$?
            if [ $RESULT -eq 0 ] ; then
                break
            fi
            echo "Transfer of '$FILE' failed. Will retry in 1 minute..."
            sleep 60
        done
        n=$(($n + 1))
        if [[ $MAX_N_FILES != 0 && $n >= $MAX_N_FILES ]] ; then
            "Max number of files was transfered. Exiting...
            exit
        fi
    done
done
