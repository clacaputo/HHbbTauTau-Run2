#!/bin/bash

if [ $# -le 1 -o $# -ge 4 ] ; then
    echo "Usage: name_file path_file_list path_output"
    exit
fi

NAME=$1
PATH_FILE_LIST=$2
PATH_OUTPUT=$3

cmsRun TreeProduction/python/treeProducer.py globalTag=START53_V7A::All includeSim=False \
  fileList=$PATH_FILE_LIST/${NAME}.txt maxEvents=-1 outputFile=$PATH_OUTPUT/${NAME}_Tree.root \
  > $PATH_OUTPUT/${NAME}_detail.log 2> $PATH_OUTPUT/${NAME}.log &
