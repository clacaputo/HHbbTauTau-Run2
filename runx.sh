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
EXE_NAME="$SCRIPT_RUN_PATH/${NAME}_${SUFFIX}"

$SCRIPT_PATH/RunTools/make.sh $SCRIPT_RUN_PATH $SUFFIX $*

RESULT=$?
if [[ $RESULT -eq 0 ]] ; then
        $EXE_NAME
	rm -f "$EXE_NAME"
fi

