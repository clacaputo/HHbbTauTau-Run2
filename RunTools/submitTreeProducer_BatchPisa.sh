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

cd $WORKING_PATH/FILE_LIST_PATH
for f in *.txt ; do
    NAME="${f%.*}"
    echo bsub -q local $RUN_SCRIPT_PATH $NAME $WORKING_PATH $FILE_LIST_PATH $OUTPUT_PATH \
                                        $GLOBAL_TAG $INCLUDE_SIM $N_EVENTS
done
