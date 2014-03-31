#!/bin/bash

if [ $# -ne 6 ] ; then
    echo "Usage: queue max_n_parallel_jobs file_list_path output_path global_tag include_sim"
    exit
fi

QUEUE=$1
MAX_N_PARALLEL_JOBS=$2
FILE_LIST_PATH=$3
OUTPUT_PATH=$4
GLOBAL_TAG=$5
INCLUDE_SIM=$6

WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
RUN_SCRIPT_PATH=$WORKING_PATH/RunTools/runTreeProducer.sh
N_EVENTS=-1

if [ ! -d "$WORKING_PATH" ] ; then
	echo "ERROR: working path '$WORKING_PATH' does not exist."
	exit
fi

if [ ! -d "$WORKING_PATH/$FILE_LIST_PATH" ] ; then
	echo "ERROR: file list path '$WORKING_PATH/$FILE_LIST_PATH' does not exist."
	exit
fi

if [ ! -d "$WORKING_PATH/$OUTPUT_PATH" ] ; then
	echo "ERROR: output path '$WORKING_PATH/$OUTPUT_PATH' does not exist."
	exit
fi

if [ ! -f "$RUN_SCRIPT_PATH" ] ; then
	echo "ERROR: script '$RUN_SCRIPT_PATH' does not exist."
	exit
fi

JOBS=$( find $WORKING_PATH/$FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" | sed "s/\.txt//" )

if [ "x$JOBS" = "x" ] ; then
	echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
	exit
fi

N_JOBS=$( echo "$JOBS" | wc -l )
echo "Following jobs will be submited:" $JOBS
echo "Total number of jobs to submit: $N_JOBS"

read -p "Submit these jobs (yes/no)? " -r REPLAY
if [ "$REPLAY" != "y" -a "$REPLAY" != "yes" -a "$REPLAY" != "Y" ] ; then 
	echo "No jobs have been submited."
	exit
fi

i=0

if [ "$QUEUE" = "local" ] ; then
    for NAME in $JOBS ; do
        bsub -q $QUEUE -J $NAME $RUN_SCRIPT_PATH $NAME $WORKING_PATH $FILE_LIST_PATH $OUTPUT_PATH \
                                            $GLOBAL_TAG $INCLUDE_SIM $N_EVENTS
    done
    echo "$N_JOBS have been submited in local"
elif [ "$QUEUE" = "fai5" -o "$QUEUE" = "fai" ] ; then
    for NAME in $JOBS ; do
        bsub -Is -q $QUEUE -J $NAME $RUN_SCRIPT_PATH $NAME $WORKING_PATH $FILE_LIST_PATH $OUTPUT_PATH \
                                                $GLOBAL_TAG $INCLUDE_SIM $N_EVENTS
        i=$(($i + 1))
        if [[ $i == $MAX_N_PARALLEL_JOBS ]] ; then
                wait
                i=0
        fi
    done
    wait
    echo "$MAX_N_PARALLEL_JOBS finished on fai"
fi
