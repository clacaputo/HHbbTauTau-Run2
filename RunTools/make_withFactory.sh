#!/bin/bash
#
#  \file make_withFactory.sh
#  \brief Compile program executable for a given C++ file, creating the main object using make_tools::Factory class.
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

SCRIPT_RUN_PATH=$1
JOB_NAME=$2
NAME=$3

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

CODE_OUT="$SCRIPT_RUN_PATH/${JOB_NAME}.cpp"
EXE_NAME="$SCRIPT_RUN_PATH/${JOB_NAME}"
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
        $( root-config --libs ) -lMathMore -lGenVector -lASImage \
        -o $EXE_NAME $CODE_OUT

RESULT=$?
rm -f "$CODE_OUT"
exit $RESULT
