#!/bin/bash

SCRIPT_RUN_PATH=$1
NAME=$2

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

CODE_OUT="$SCRIPT_RUN_PATH/${NAME}.cpp"
EXE_NAME="$SCRIPT_RUN_PATH/${NAME}"
rm -f "$EXE_NAME"

printf \
"#include \"${SOURCE}\"
#include <TROOT.h>
#include <iostream>
#include <memory>
int main(int argc, char *argv[])
{
        try {
                gROOT->ProcessLine(\"#include <vector>\");
                std::unique_ptr<$NAME> a( make_tools::Factory<$NAME>::Make(argc, argv) );
                a->Run();
        }
        catch(std::exception& e) {
                std::cerr << \"ERROR: \" << e.what() << std::endl;
                return 1;
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
