#!/bin/bash

if [ $# -ne 2 ] ; then
    echo "Usage: dir_name file_name_prefix"
    echo "Output: paths of files that has the same job id."
    exit
fi

DIR_NAME=$1
PREFIX=$2

ls $DIR_NAME | \
    sed "s/\(${PREFIX}_\)\([0-9]*\)\(.*\)/\2/" | sort | uniq -d | \
    sed "s/\(.*\)/${PREFIX}_\1_\*/" | \
    xargs -n 1 find $DIR_NAME -name
