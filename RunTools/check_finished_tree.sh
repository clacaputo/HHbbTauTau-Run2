#!/bin/bash

if [ $# -ne 1 ] ; then
    echo "Usage: file_list_path "
    exit
fi

FILE_LIST_PATH=$1


WORKING_PATH=$PWD

if [ ! -d "$WORKING_PATH" ] ; then
        echo "ERROR: working path '$WORKING_PATH' does not exist."
        exit
fi

if [ ! -d "$WORKING_PATH/$FILE_LIST_PATH" ] ; then
        echo "ERROR: file list path '$WORKING_PATH/$FILE_LIST_PATH' does not exist."
        exit
fi


JOBS=$( find $WORKING_PATH/$FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" )

if [ "x$JOBS" = "x" ] ; then
        echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
        exit
fi

N_JOBS=$( echo "$JOBS" | wc -l )
echo "Following jobs will be submited:" $JOBS
echo "Total number of jobs to submit: $N_JOBS"
