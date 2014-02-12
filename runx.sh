#!/bin/bash

function get_arg_value {
	if [ "${ARGS[$1]:0:1}" == "@" ] ; then
		arg_value="${ARGS[$1]:1}"
	else
		arg_value="\"${ARGS[$1]}\""
	fi
	
}

NAME=$1
ARGS=($*)
SCRIPT_PATH="$CMSSW_BASE/src/HHbbTauTau"
SCRIPT_RUN_PATH="$SCRIPT_PATH/run"
mkdir -p $SCRIPT_RUN_PATH

SOURCE=$( find $SCRIPT_PATH -name "${NAME}.C" )
n=${#ARGS[@]}
arg_list=""
if (( n > 1 )) ; then
	get_arg_value 1
	arg_list="$arg_value"
	for (( i = 2; i < n; i++ )) ; do
		get_arg_value i
		arg_list="$arg_list, $arg_value"
	done
fi

CODE_OUT="$SCRIPT_RUN_PATH/${NAME}.cpp"
EXE_NAME="$SCRIPT_RUN_PATH/$NAME"
rm -f $EXE_NAME

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

g++ \
	-Isrc/ -I$CMSSW_RELEASE_BASE/src -I$ROOT_INCLUDE_PATH -I$BOOST_INCLUDE_PATH \
	`root-config --libs` \
	-o $EXE_NAME $CODE_OUT

$EXE_NAME
