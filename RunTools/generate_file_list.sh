#!/bin/bash

if [ $# -le 1 -o $# -ge 4 ] ; then
    echo "Usage: input_path output_file_name [reference_path]"
    exit
fi

INPUT_PATH=$1
OUT_FILE=$(echo $(cd $(dirname $2); pwd)/$(basename $2))
REFERENCE_PATH=$3
PREFIX=""

if [ $# -ge 3 ] ; then
    cd "$REFERENCE_PATH"
    PREFIX="/"
    if [ "${INPUT_PATH:0:1}" == "/" ] ; then
        INPUT_PATH=${INPUT_PATH:1}
    fi
fi

find  $INPUT_PATH -name "*.root" -printf "${PREFIX}%p\n" > $OUT_FILE
