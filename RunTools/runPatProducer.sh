#!/bin/bash
#
#  \file runPatProducer.sh
#  \brief Run PAT-production jobs with environment setup and output redirection.
#  \author Konstantin Androsov (University of Siena, INFN Pisa)
#  \author Maria Teresa Grippo (University of Siena, INFN Pisa)
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

if [ $# -ne 9 ] ; then
    echo "Usage: job_name working_path file_list_path output_path global_tag is_mc include_sim n_events prefix"
    exit
fi

NAME=$1
WORKING_PATH=$2
FILE_LIST_PATH=$3
OUTPUT_PATH=$4
GLOBAL_TAG=$5
IS_MC=$6
INCLUDE_SIM=$7
N_EVENTS=$8
PREFIX=$9

KEEP_PAT=False
RUN_TREE=True

cd $WORKING_PATH

echo "$NAME $( date )" >> $OUTPUT_PATH/job_start.log

source cmsenv.sh
eval $( scramv1 runtime -sh )

echo "$NAME $( date )" >> $OUTPUT_PATH/job_cmsRun_start.log

if [ $PREFIX = "none" ] ; then
    cmsRun PatProduction/python/patTuple.py globalTag=$GLOBAL_TAG isMC=$IS_MC includeSim=$INCLUDE_SIM \
           fileList=$FILE_LIST_PATH/${NAME}.txt maxEvents=$N_EVENTS outputFile=$OUTPUT_PATH/${NAME}_Pat.root \
           keepPat=$KEEP_PAT runTree=$RUN_TREE treeOutput=$OUTPUT_PATH/${NAME}_Tree.root \
           > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
    RESULT=$?
else
    cmsRun PatProduction/python/patTuple.py globalTag=$GLOBAL_TAG isMC=$IS_MC includeSim=$INCLUDE_SIM \
           fileList=$FILE_LIST_PATH/${NAME}.txt maxEvents=$N_EVENTS outputFile=$OUTPUT_PATH/${NAME}_Pat.root \
           keepPat=$KEEP_PAT runTree=$RUN_TREE treeOutput=$OUTPUT_PATH/${NAME}_Tree.root \
           fileNamePrefix=$PREFIX > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
    RESULT=$?
fi


echo "$RESULT $NAME $( date )" >> $OUTPUT_PATH/job_result.log

exit $RESULT
