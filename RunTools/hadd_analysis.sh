#!/bin/bash

if [ $# -ne 2 ] ; then
    echo "Usage: analysis_path output_path"
    exit
fi

ANALYSIS_PATH=$1
