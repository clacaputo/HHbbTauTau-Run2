#!/bin/bash
#
#  \file batch.sh
#  \brief Definition of a function that submit jobs on batch.
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

function submit_batch {

    if [ $# -lt 5 ] ; then
        echo "Usage: queue storage job_name output_path script [script_args]"
        exit 1
    fi

    QUEUE=$1
    STORAGE=$2
    JOB_NAME=$3
    OUTPUT_PATH=$4
    RUN_SCRIPT_PATH=$5

    if [ ! -f "$RUN_SCRIPT_PATH" ] ; then
        echo "ERROR: script '$RUN_SCRIPT_PATH' does not exist."
        exit 1
    fi

    if [ "$QUEUE" = "cms" -a "$STORAGE" = "Pisa" ] ; then
        bsub -q $QUEUE -m austroh -E /usr/local/lsf/work/infn-pisa/scripts/testq-preexec-cms.bash \
             -J $JOB_NAME $RUN_SCRIPT_PATH ${@:6}
    elif [ "$QUEUE" = "local" -a "$STORAGE" = "Pisa" ] ; then
        bsub -q $QUEUE -m austroh -E /usr/local/lsf/work/infn-pisa/scripts/testq-preexec-cms.bash \
             -J $JOB_NAME $RUN_SCRIPT_PATH ${@:6}
    elif [ "$QUEUE" = "local" -a "$STORAGE" = "Bari" ] ; then
        echo "$RUN_SCRIPT_PATH ${@:6}" | qsub -q $QUEUE -N $JOB_NAME -o $OUTPUT_PATH -e $OUTPUT_PATH -
    elif [ "$QUEUE" = "fai5" -a "$STORAGE" = "Pisa" ] ; then
        bsub -Is -q $QUEUE -J $JOB_NAME $RUN_SCRIPT_PATH ${@:6} &
    elif [ "$QUEUE" = "fai" ] ; then
        bsub -Is -q $QUEUE -R "select[defined(fai)]" -J $JOB_NAME $RUN_SCRIPT_PATH ${@:6} &
    elif [ $STORAGE = "Local" ] ; then
        $RUN_SCRIPT_PATH ${@:6} &
    else
        echo "unknow queue"
        exit 1
    fi
}

