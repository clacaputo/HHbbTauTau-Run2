#!/bin/bash

if [ $# -ne 2 ] ; then
    echo "Usage: trees_path even_id"
    exit
fi

TREES_PATH=$1
EVENT_ID=$2

SCAN_FILE="RunTools/source/Scan.C"

find $TREES_PATH -type f -exec root -b -l -q $SCAN_FILE\($EVENT_ID,\"\{\}\"\) \; > /dev/null
