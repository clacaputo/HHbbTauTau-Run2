#!/bin/bash

if [ $# -ne 7 ] ; then
    echo "Usage: job_name working_path file_list_path output_path global_tag include_sim n_events"
    exit
fi

NAME=$1
WORKING_PATH=$2
FILE_LIST_PATH=$3
OUTPUT_PATH=$4
GLOBAL_TAG=$5
INCLUDE_SIM=$6
N_EVENTS=$7

cd $WORKING_PATH
source cmsenv.sh
eval $( scramv1 runtime -sh )

echo "$NAME $( date )" >> $OUTPUT_PATH/job_start.log

cmsRun TreeProduction/python/treeProducer.py globalTag=$GLOBAL_TAG includeSim=$INCLUDE_SIM \
       fileList=$FILE_LIST_PATH/${NAME}.txt maxEvents=$N_EVENTS outputFile=$OUTPUT_PATH/${NAME}_Tree.root \
       > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
RESULT=$?

echo "$RESULT $NAME $( date )" >> $OUTPUT_PATH/job_result.log

exit $RESULT
