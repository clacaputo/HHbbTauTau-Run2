#!/bin/bash
#
#  \file runAnalysis.sh
#  \brief Run analysis jobs with environment setup and output redirection.
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

if [ $# -lt 5 ] ; then
    echo "Usage: job_name working_path output_path exe_name set_cmsenv [other_args_for_analyzer]"
    exit
fi

NAME=$1
WORKING_PATH=$2
OUTPUT_PATH=$3
EXE_NAME=$4
SET_CMSENV=$5

cd $WORKING_PATH

echo "$NAME $( date )" >> $OUTPUT_PATH/job_start.log

if [ $SET_CMSENV = "yes" ] ; then
    source cmsenv.sh
    eval $( scramv1 runtime -sh )
    echo "$NAME $( date )" >> $OUTPUT_PATH/job_cmsRun_start.log
fi

$EXE_NAME "${@:6}" > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
RESULT=$?

echo "$RESULT $NAME $( date )" >> $OUTPUT_PATH/job_result.log

exit $RESULT
