#!/bin/bash

function get_arg_value {
        if [ "${ARGS[$1]:0:1}" == "@" ] ; then
                arg_value="${ARGS[$1]:1}"
        else
                arg_value="\"${ARGS[$1]}\""
        fi

}


SCRIPT_RUN_PATH=$1
JOB_NAME=$2
NAME=$3
ARGS=($*)

if [ "$CMSSW_BASE/" = "/" ] ; then
    SCRIPT_PATH="."
else
    SCRIPT_PATH="$CMSSW_BASE/src/HHbbTauTau"
fi

if [ ! -d "$SCRIPT_RUN_PATH" ] ; then
    echo "ERROR: script run path $SCRIPT_RUN_PATH doesn't exist"
    exit 1
fi

SOURCE=$( find $SCRIPT_PATH -name "${NAME}.C" )
n=${#ARGS[@]}
arg_list=""
if (( n > 3 )) ; then
        get_arg_value 3
        arg_list="$arg_value"
        for (( i = 4; i < n; i++ )) ; do
                get_arg_value i
                arg_list="$arg_list, $arg_value"
        done
fi

CODE_OUT="$SCRIPT_RUN_PATH/${JOB_NAME}.cpp"
EXE_NAME="$SCRIPT_RUN_PATH/${JOB_NAME}"
rm -f "$EXE_NAME"

printf \
"#include \"${SOURCE}\"
#include <TROOT.h>
#include <iostream>
int main()
{
        try {
                gROOT->ProcessLine(\"#include <vector>\");
                $NAME a( $arg_list );
                a.Run();
        }
        catch(std::exception& e) {
                std::cerr << \"ERROR: \" << e.what() << std::endl;
        }
        return 0;
}
" > $CODE_OUT

g++ -std=c++0x -Wall -O3 \
        -I. -I$CMSSW_BASE/src -I$CMSSW_RELEASE_BASE/src -I$ROOT_INCLUDE_PATH -I$BOOST_INCLUDE_PATH \
        $( root-config --libs ) -lMathMore -lGenVector \
        -o $EXE_NAME $CODE_OUT

RESULT=$?
rm -f "$CODE_OUT"
exit $RESULT
