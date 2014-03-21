#!/bin/bash

#if [ $# -le 1 -o $# -ge 6 ] ; then
#    echo "Usage: name_file path_file_list path_output n_events"
#    exit
#fi

#NAME=$1
#PATH_FILE_LIST=$2
#PATH_OUTPUT=$3
#N_EVENTS=$4

NAME=DYtautau_xaa 
PATH_FILE_LIST=TreeProduction/dataset/DYtautau_Pat 
PATH_OUTPUT=data/DYtautau
N_EVENTS=1

ANALYSIS_DIR=/gpfs/ddn/cms/user/grippo/HHbbtautau/CMSSW_5_3_14_patch2/src/HHbbTauTau
cd $ANALYSIS_DIR
. cmsenv.sh
eval $( scramv1 runtime -sh )

cmsRun TreeProduction/python/treeProducer.py globalTag=START53_V7A::All includeSim=False \
  fileList=$PATH_FILE_LIST/${NAME}.txt maxEvents=$N_EVENTS outputFile=$PATH_OUTPUT/${NAME}_Tree.root \
  > $PATH_OUTPUT/${NAME}_detail.log 2> $PATH_OUTPUT/${NAME}.log
