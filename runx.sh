#!/bin/bash

NAME=$1

if [ "$CMSSW_BASE/" = "/" ] ; then
    SCRIPT_PATH="."
else
    SCRIPT_PATH="$CMSSW_BASE/src/HHbbTauTau"
fi

SCRIPT_RUN_PATH="$SCRIPT_PATH/run"
mkdir -p $SCRIPT_RUN_PATH


SUFFIX=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
JOB_NAME=${NAME}_${SUFFIX}
EXE_NAME="$SCRIPT_RUN_PATH/$JOB_NAME"

$SCRIPT_PATH/RunTools/make.sh $SCRIPT_RUN_PATH $JOB_NAME $*

RESULT=$?
if [[ $RESULT -eq 0 ]] ; then
        $EXE_NAME
	rm -f "$EXE_NAME"
fi

