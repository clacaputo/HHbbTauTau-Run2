#!/bin/bash

if [ $# -ne 1 ] ; then
    echo "Usage: reference_input_path output_path "
    exit
fi

REFERENCE_INPUT_PATH=$1

if [ ! -d "$REFERENCE_INPUT_PATH" ] ; then
    echo "ERROR: working path '$REFERENCE_INPUT_PATH' does not exist."
    exit
fi

OUTPUT_PATH=$2

if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: working path '$OUTPUT_PATH' does not exist."
    exit
fi

DATASET_LIST=$( ls $REFERENCE_INPUT_PATH )

ANALYZER_PATH="HHbbMuTaujet_analyzer HHbbETaujet_analyzer HHbb2Taujet_analyzer"
WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
RUN_SCRIPT_PATH=$WORKING_PATH/RunTools/submitAnalysis_Batch.sh

if [ ! -d "$WORKING_PATH" ] ; then
    echo "ERROR: working path '$WORKING_PATH' does not exist."
    exit
fi

for ANA_FOLDER in $ANALYZER_PATH ; do
    mkdir -p $OUTPUT_PATH/$ANA_FOLDER

    if [ -d "$OUTPUT_PATH/$ANA_FOLDER" ] ; then
        echo "ERROR: output dir '$OUTPUT_PATH/$ANA_FOLDER' already exists."
        exit
    fi
    for DATASET_DIR in $DATASET_LIST ; do
        INPUT_PATH=$REFERENCE_INPUT_PATH/$DATASET_DIR
        mkdir -p $OUTPUT_PATH/$ANA_FOLDER/$DATASET_DIR

        if [ -d "$OUTPUT_PATH/$ANA_FOLDER/$DATASET_DIR" ] ; then
            echo "ERROR: output dir '$OUTPUT_PATH/$ANA_FOLDER/$DATASET_DIR' already exists."
            exit
        fi
        $RUN_SCRIPT_PATH local Bari 0 $INPUT_PATH $OUTPUT_PATH/$ANA_FOLDER/$DATASET_DIR \
                         $ANA_FOLDER @false Analysis/data/reweight_PU.root
    done
done


