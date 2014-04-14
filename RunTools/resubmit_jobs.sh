#!/bin/bash

if [ $# -ne 5 ] ; then
    echo "Usage: dataset_name output_path global_tag include_sim prefix"
    exit
fi

DATASET=$1
OUTPUT_PATH=$2

if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: output path '$OUTPUT_PATH' does not exist."
    exit
fi
OUTPUT_PATH=$( cd "$OUTPUT_PATH" ; pwd )


FILE_LIST_PATH="TreeProduction/dataset/${DATASET}_Pat"
if [ ! -d $FILE_LIST_PATH ] ; then
    echo "ERROR: file list path '$FILE_LIST_PATH' does not exists."
    exit
fi
JOBS=$( find $FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" | sed "s/\.txt//" | sort )
if [ "x$JOBS" = "x" ] ; then
        echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
        exit
fi

NEW_FILE_LIST_PATH="TreeProduction/dataset/retry/${DATASET}"
n=1
while [-d $NEW_FILE_LIST_PATH ] ; do
    n=(($n + 1))
    NEW_FILE_LIST_PATH="TreeProduction/dataset/retry/${DATASET}_${n}"
done
echo "Resubmit number $n"
mkdir -p $NEW_FILE_LIST_PATH

FILE_JOB_RESULT="$OUTPUT_PATH/job_result.log"
SUCCESSFULL_JOBS=""
N_SUCCESSFULL_JOBS=0
if [ -f $FILE_JOB_RESULT ] ; then
    SUCCESSFULL_JOBS=$( cat $FILE_JOB_RESULT | sed -n "s/\(^0 \)\([^ ]*\)\(.*\)/\2/p" | sort )
    N_SUCCESSFULL_JOBS=$( echo "$SUCCESSFULL_JOBS" | wc -l )
    echo "Number of successfull jobs: $N_SUCCESSFULL_JOBS"
else
    echo "job_results.log not found. Considering that there are no finished jobs yet."
fi

JOB_STAT_RESULT=$( qstat -u $(whoami) | grep $DATASET )
JOBS_IN_QUEUE=$( echo "$JOB_STAT_RESULT" | awk '{print $4}' )
N_RUNNING_JOBS=$( echo "$JOB_STAT_RESULT" | grep " R " | wc -l )
N_PENDING_JOBS=$( echo "$JOB_STAT_RESULT" | grep " Q " | wc -l )

echo "Number of running jobs: $N_RUNNING_JOBS"
echo "Number of pending jobs: $N_PENDING_JOBS"

echo -e "$JOBS"\\n"$SUCCESSFULL_JOBS"\\n"$JOBS_IN_QUEUE" | sort | uniq -u \
    | xargs -n 1 printf "$FILE_LIST_PATH/%b.txt $NEW_FILE_LIST_PATH/\n" | xargs -n 2 cp

RunTools/submitTreeProducer_Batch.sh local Bari 0 $NEW_FILE_LIST_PATH $OUTPUT_PATH $GLOBAL_TAG $INCLUDE_SIM $PREFIX
