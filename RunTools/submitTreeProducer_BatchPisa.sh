#!/bin/bash

if [ $# -ne 4 ] ; then
    echo "Usage: file_list_path output_path global_tag include_sim"
    exit
fi

FILE_LIST_PATH=$1
OUTPUT_PATH=$2
GLOBAL_TAG=$3
INCLUDE_SIM=$4

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

for NAME in $JOBS ; do
    bsub -q local -J $NAME $RUN_SCRIPT_PATH $NAME $WORKING_PATH $FILE_LIST_PATH $OUTPUT_PATH \
                                            $GLOBAL_TAG $INCLUDE_SIM $N_EVENTS
done

echo "$N_JOBS have been submited."
