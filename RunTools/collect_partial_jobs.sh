#!/bin/bash

if [ $# -le 2 -o $# -ge 5 ] ; then
    echo "Usage: input_path output_file_name n_jobs [reference_path]"
    exit
fi

INPUT_PATH=$1
OUT_DIR=$(echo $(cd $(dirname $2); pwd))
OUT_FILE_NAME=$(basename $2)
OUT_FILE=$OUT_DIR/$OUT_FILE_NAME
TMP_DIR=$(mktemp -d)
TMP_OUT_FILE=$TMP_DIR/$OUT_FILE_NAME
N_JOBS=$3
REFERENCE_PATH=$4
PREFIX=""

if [ $# -ge 3 ] ; then
    cd "$REFERENCE_PATH"
    PREFIX="/"
    if [ "${INPUT_PATH:0:1}" == "/" ] ; then
        INPUT_PATH=${INPUT_PATH:1}
    fi
fi

FILE_JOB_RESULT="$INPUT_PATH/job_result.log"
if [ -f $FILE_JOB_RESULT ] ; then
    echo "ERROR: job_result.log not found."
    exit
fi

SUCCESSFULL_JOBS=$( cat $FILE_JOB_RESULT | sed -n "s/\(^0 \)\([^ ]*\)\(.*\)/\2/p" | sort )
N_SUCCESSFULL_JOBS=$( echo "$SUCCESSFULL_JOBS" | sed '/^\s*$/d' | wc -l )
if [ $N_SUCCESSFULL_JOBS -eq 0 ] ; then
    echo "There are no successfull jobs found. Exiting..."
    exit
fi
echo "Number of successfull jobs: $N_SUCCESSFULL_JOBS"

echo "$SUCCESSFULL_JOBS" | xargs -n 1 printf "${PREFIX}$INPUT_PATH/%b_Tree.root" > $TMP_OUT_FILE

mkdir -p $OUTPUT_DIR

cd $TMP_DIR
N_FILES=$( cat $TMP_OUT_FILE | wc -l )
N_SPLIT=$(( (N_FILES + N_JOBS - 1) / N_JOBS ))
echo "Total number of files: $N_FILES"
echo "Total number of jobs: $N_JOBS"
echo "Number of files per job: $N_SPLIT"

if [ $N_JOBS -eq 1 ] ; then
    mv "$TMP_OUT_FILE" "${OUT_FILE_NAME}.txt"
else
    split -l $N_SPLIT $TMP_OUT_FILE
    find . -maxdepth 1 -type f -name "x*" -printf "%f ${OUT_FILE_NAME}_%f.txt\n" | xargs -n 2 mv
fi

mv *.txt $OUT_DIR
rm -rf $TMP_DIR
