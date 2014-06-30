#!/bin/bash
#
#  \file submitAllAnalysis.sh
#  \brief Submit all available analysis datasets on batch.
#  \author Konstantin Androsov (Siena University, INFN Pisa)
#  \author Maria Teresa Grippo (Siena University, INFN Pisa)
#
#  Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
#                 Maria Teresa Grippo <grippomariateresa@gmail.com>
#
#  This file is part of X->HH->bbTauTau.
#
#  X->HH->bbTauTau is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  X->HH->bbTauTau is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.

if [ $# -ne 7 ] ; then
    echo "Usage: analyzer_cfg dataset_cfg output_path queue storage n_parallel_jobs use_multi_ana"
    exit
fi

ANA_CFG_FILE=$1
CFG_FILE=$2

if [ ! -f "$CFG_FILE" ] ; then
    echo "ERROR: config file '$CFG_FILE' does not exist."
    exit
fi

OUTPUT_PATH=$3

if [ -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: working path '$OUTPUT_PATH' already exists."
    exit
fi
OUTPUT_PATH=$( cd "$OUTPUT_PATH" ; pwd )

QUEUE=$4
STORAGE=$5
MAX_N_PARALLEL_JOBS=$6
USE_MULTI_ANA=$7

WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
MAKE_PATH=$WORKING_PATH/RunTools/make_withFactory.sh

REFERENCE_INPUT_PATH="Analysis/dataset"
REFERENCE_CFG_PATH="Analysis/config"
DATASET_ARRAY=( $( cat $CFG_FILE | awk '{ print $1 }' ) )
CFG_ARRAY=( $( cat $CFG_FILE | awk '{ print $2 }' ) )

N_DATASET=${#DATASET_ARRAY[@]}

ANALYZER_PATH=( $( cat $ANA_CFG_FILE ) )

if [ $USE_MULTI_ANA = "yes" ] ; then
    MULTI_ANA_NAME="H_BaselineSync"

    $MAKE_PATH $OUTPUT_PATH $MULTI_ANA_NAME $MULTI_ANA_NAME
    EXE_NAME=$OUTPUT_PATH/$MULTI_ANA_NAME/$MULTI_ANA_NAME
    echo "Executable file $EXE_NAME is compiled."

    for (( i=0; i<$N_DATASET; i++ )) ; do
        INPUT_PATH=$REFERENCE_INPUT_PATH/${DATASET_ARRAY[i]}

        if [ ! -d "$INPUT_PATH" ] ; then
            echo "ERROR: dataset '${DATASET_ARRAY[i]}' does not exist."
            exit
        fi

        ANA_OUTPUT_PATH_1=$OUTPUT_PATH/${ANALYZER_PATH[0]}/${DATASET_ARRAY[i]}
        mkdir -p $ANA_OUTPUT_PATH_1

        ANA_OUTPUT_PATH_2=$OUTPUT_PATH/${ANALYZER_PATH[1]}/${DATASET_ARRAY[i]}
        mkdir -p $ANA_OUTPUT_PATH_2

        ANA_OUTPUT_PATH_3=$OUTPUT_PATH/${ANALYZER_PATH[2]}/${DATASET_ARRAY[i]}
        mkdir -p $ANA_OUTPUT_PATH_3

        yes | ./RunTools/submitMultiAnalysis_Batch.sh $QUEUE $STORAGE $MAX_N_PARALLEL_JOBS $INPUT_PATH \
            $ANA_OUTPUT_PATH_1 $ANA_OUTPUT_PATH_2 $ANA_OUTPUT_PATH_3 \
            $ANA_FOLDER $REFERENCE_CFG_PATH/${CFG_ARRAY[i]} $EXE_NAME
    done

else

    for ANA_FOLDER in $ANALYZER_PATH ; do
        mkdir -p $OUTPUT_PATH/$ANA_FOLDER

        $MAKE_PATH $OUTPUT_PATH/$ANA_FOLDER $ANA_FOLDER $ANA_FOLDER
        EXE_NAME=$OUTPUT_PATH/$ANA_FOLDER/$ANA_FOLDER
        echo "Executable file $EXE_NAME is compiled."

        for (( i=0; i<$N_DATASET; i++ )) ; do
            INPUT_PATH=$REFERENCE_INPUT_PATH/${DATASET_ARRAY[i]}

            if [ ! -d "$INPUT_PATH" ] ; then
                echo "ERROR: dataset '${DATASET_ARRAY[i]}' does not exist."
                exit
            fi

            ANA_OUTPUT_PATH=$OUTPUT_PATH/$ANA_FOLDER/${DATASET_ARRAY[i]}
            mkdir -p $ANA_OUTPUT_PATH

            echo $ANA_FOLDER ${DATASET_ARRAY[i]} ${CFG_ARRAY[i]}

            yes | ./RunTools/submitAnalysis_Batch.sh $QUEUE $STORAGE $MAX_N_PARALLEL_JOBS $INPUT_PATH $ANA_OUTPUT_PATH \
                $ANA_FOLDER $REFERENCE_CFG_PATH/${CFG_ARRAY[i]} $EXE_NAME

        done

    done

fi

