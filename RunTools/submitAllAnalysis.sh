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

if [ $# -ne 2 ] ; then
    echo "Usage: reference_input_path output_path "
    exit
fi

REFERENCE_INPUT_PATH=$1

if [ ! -d "$REFERENCE_INPUT_PATH" ] ; then
    echo "ERROR: working path '$REFERENCE_INPUT_PATH' does not exist."
    exit
fi

OUTPUT_PATH=$2

if [ -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: working path '$OUTPUT_PATH' already exists."
    exit
fi

DATASET_LIST=$( ls $REFERENCE_INPUT_PATH )

ANALYZER_PATH="HHbbMuTaujet_analyzer HHbbETaujet_analyzer HHbb2Taujet_analyzer"

for ANA_FOLDER in $ANALYZER_PATH ; do
    mkdir -p $OUTPUT_PATH/$ANA_FOLDER

    for DATASET_DIR in $DATASET_LIST ; do
        INPUT_PATH=$REFERENCE_INPUT_PATH/$DATASET_DIR
        mkdir -p $OUTPUT_PATH/$ANA_FOLDER/$DATASET_DIR

        RunTools/submitAnalysis_Batch.sh local Bari 0 $INPUT_PATH $OUTPUT_PATH/$ANA_FOLDER/$DATASET_DIR \
                         $ANA_FOLDER @false Analysis/data/reweight_PU.root
    done
done


