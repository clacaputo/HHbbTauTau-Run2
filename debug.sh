#!/bin/bash
#
#  \file debug.sh
#  \brief Compile and run C++ file with a specified name, passing to it command line arguments in debug mode.
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

NAME=$1

if [ "$CMSSW_BASE/" = "/" ] ; then
    SCRIPT_PATH="."
else
    SCRIPT_PATH="$CMSSW_BASE/src/HHbbTauTau"
fi

SCRIPT_RUN_PATH="$SCRIPT_PATH/run"
mkdir -p $SCRIPT_RUN_PATH


LC_CTYPE=C SUFFIX=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
JOB_NAME=${NAME}_${SUFFIX}
EXE_NAME="$SCRIPT_RUN_PATH/$JOB_NAME"

$SCRIPT_PATH/RunTools/make.sh $SCRIPT_RUN_PATH $JOB_NAME -g $*

RESULT=$?
if [[ $RESULT -eq 0 ]] ; then
        valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=debug.log \
                 --suppressions=RunTools/config/known_memory_leaks.supp \
                 --suppressions=RunTools/config/known_issues.supp --gen-suppressions=all $EXE_NAME
	rm -f "$EXE_NAME"
fi

